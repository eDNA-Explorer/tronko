/*
 *	getclade.c
 */
#include "getclade.h"
#include <string.h>

/* ---- leaf name -> index lookup -------------------------------------------
   The leaf branch below used to strcmp every entry of m->names for every leaf
   in the tree, and did not stop at a match -- O(numspec^2 * namelen) per parse,
   run once per partition at every level of the recursion. For the largest CO1
   cluster (182,587 leaves) that is 3.3e10 strcmp calls per parse.

   Open-addressed table, built once per parse. Insert is last-wins so a
   duplicate name resolves to the same index the original scan would have left
   in `tip` (it kept overwriting rather than breaking).

   __thread because lane worker pthreads parse concurrently -- same reason
   global.h:39 declares root/tip/comma as __thread. */
static __thread const char **nm_key;
static __thread int         *nm_val;
static __thread int          nm_cap;   /* power of two; 0 = not built */

static unsigned long nm_hash(const char *s){
	unsigned long h = 1469598103934665603UL;      /* FNV-1a */
	while (*s){ h ^= (unsigned char)*s++; h *= 1099511628211UL; }
	return h;
}
static void nm_build(struct masterArr *m){
	int i;
	nm_cap = 1;
	while (nm_cap < m->numspec * 2) nm_cap <<= 1;
	nm_key = calloc(nm_cap, sizeof(char *));
	nm_val = malloc(nm_cap * sizeof(int));
	if (!nm_key || !nm_val){            /* fall back to the linear scan */
		free(nm_key); free(nm_val);
		nm_key = NULL; nm_val = NULL; nm_cap = 0;
		return;
	}
	for (i = 0; i < m->numspec; i++){
		unsigned long h = nm_hash(m->names[i]) & (unsigned long)(nm_cap - 1);
		while (nm_key[h] && strcmp(nm_key[h], m->names[i]) != 0)
			h = (h + 1) & (unsigned long)(nm_cap - 1);
		nm_key[h] = m->names[i];
		nm_val[h] = i;
	}
}
static void nm_free(void){
	free(nm_key); free(nm_val);
	nm_key = NULL; nm_val = NULL; nm_cap = 0;
}
/* -1 = not found, -2 = table unavailable (caller must scan) */
static int nm_lookup(const char *s){
	unsigned long h;
	if (!nm_cap) return -2;
	h = nm_hash(s) & (unsigned long)(nm_cap - 1);
	while (nm_key[h]){
		if (strcmp(nm_key[h], s) == 0) return nm_val[h];
		h = (h + 1) & (unsigned long)(nm_cap - 1);
	}
	return -1;
}

/* O(n) in-memory Newick parser — replaces O(n^2) FILE-based getcladeArr.
   Reads binary Newick string from memory buffer using pointer arithmetic.
   No fgetc/fsetpos overhead. */
int getcladeArr_mem(const char **pp, struct masterArr *m, int max_nodename){
	const char *p = *pp;
	/* Skip whitespace */
	while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') p++;

	if (*p == '('){
		/* Internal node */
		p++; /* skip '(' */
		int n1 = getcladeArr_mem(&p, m, max_nodename);
		/* skip comma separator between children */
		while (*p == ',' || *p == ' ' || *p == '\t' || *p == '\n' || *p == '\r'){
			if (*p == ',') comma++;
			p++;
		}
		/* Node number must be captured right after separator comma,
		   matching the original getnodenumbArr formula */
		int n3 = comma;
		int n2 = getcladeArr_mem(&p, m, max_nodename);
		/* skip any trailing commas/whitespace before ')' */
		while (*p == ',' || *p == ' ' || *p == '\t' || *p == '\n' || *p == '\r'){
			if (*p == ',') comma++;
			p++;
		}
		if (*p == ')') p++;
		/* Skip optional support value (e.g. "0.996" in ")0.996:0.068") */
		if (*p != ':' && *p != ',' && *p != ')' && *p != ';' && *p != '(' && *p != '\0'){
			char *endptr;
			strtod(p, (char **)&endptr);
			if (endptr != p) p = endptr;
		}
		/* Read optional branch length */
		if (*p == ':'){
			p++;
			char *endptr;
			m->tree[n3-1].bl = strtod(p, (char **)&endptr);
			p = endptr;
		}
		/* Link nodes */
		if (n1 != -1 && n2 != -1 && n3 != -1){
			linknodesArr(n1-1, n2-1, n3-1, m);
		}
		*pp = p;
		return n3;
	} else if (*p != ')' && *p != ',' && *p != ';' && *p != '\0'){
		/* Leaf node */
		char specname[max_nodename+1];
		int i = 0;
		while (*p != ':' && *p != ')' && *p != ',' && *p != ';' && *p != '\0' && i < max_nodename){
			char ch = *p;
			if (isalpha(ch) || isdigit(ch) || ch=='.' || ch=='_' || ch=='/' || ch=='-' || ch=='|' || ch=='?' || ch=='*' || ch=='&' || ch=='+' || ch=='#' || ch=='\'' || ch=='=' || ch==',' || ch==']' || ch=='[' || ch=='>' || ch=='<'){
				specname[i++] = ch;
			}
			p++;
		}
		specname[i] = '\0';
		/* Find matching species in names array (hash lookup -- see nm_build) */
		int found = 0;
		int idx = nm_lookup(specname);
		if (idx == -2){                       /* table unavailable: original scan */
			int tmp1;
			for (tmp1 = 0; tmp1 < m->numspec; tmp1++){
				if (!strcmp(m->names[tmp1], specname)){
					tip = tmp1 + 1;
					found = 1;
				}
			}
		}else if (idx >= 0){
			tip = idx + 1;
			found = 1;
		}
		if (!found) fprintf(stderr, "WARNING: leaf '%s' not found in MSA\n", specname);
		int nodeIdx = tip + m->numspec - 2;
		strcpy(m->tree[nodeIdx].name, specname);
		m->tree[nodeIdx].up[0] = -1;
		m->tree[nodeIdx].up[1] = -1;
		/* Read branch length */
		if (*p == ':'){
			p++;
			char *endptr;
			m->tree[nodeIdx].bl = strtod(p, (char **)&endptr);
			p = endptr;
		}
		*pp = p;
		return tip + (m->numspec - 1);
	}
	*pp = p;
	return -1;
}

/* Wrapper: read file into memory, parse with O(n) parser */
int getcladeArr_fast(const char *filepath, struct masterArr *m, int max_nodename){
	FILE *fp = fopen(filepath, "r");
	if (!fp){ fprintf(stderr, "Cannot open tree file: %s\n", filepath); exit(-1); }
	fseek(fp, 0, SEEK_END);
	long fsize = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	char *buf = malloc(fsize + 1);
	if (!buf){ fclose(fp); fprintf(stderr, "Out of memory reading tree\n"); exit(-1); }
	fread(buf, 1, fsize, fp);
	buf[fsize] = '\0';
	fclose(fp);

	comma = 0;
	tip = 0;
	const char *p = buf;
	nm_build(m);
	int root = getcladeArr_mem(&p, m, max_nodename);
	nm_free();
	free(buf);
	return root;
}

/* Original O(n^2) FILE-based parser — kept for compatibility */
int getcladeArr(FILE *tre, struct masterArr *m, int max_nodename){
	int n1, n2, n3;
	char ch;
	n1=-1;
	n2=-1;
	n3=-1;
	do{
		if (specsearchArr(tre,m,max_nodename)==1){
			return tip+(m->numspec-1);
		}
		ch = fgetc(tre);
		if (ch==','){
			comma++;
		}
		if (ch==')'){
			if ((ch=fgetc(tre))!=':'){
				ungetc(ch,tre);
			}else{
				do{
					ch=(fgetc(tre));
				}while((ch=='\n')||(ch==' '));
				ungetc(ch,tre);
				fscanf(tre,"%lf",&m->tree[n3-1].bl);
			}
			return n3;
		}
		if (ch=='('){
			n3=getnodenumbArr(tre);
			n1=getcladeArr(tre,m,max_nodename);
			n2=getcladeArr(tre,m,max_nodename);
			if ( n1!=-1 && n2!=-1 && n3!=-1){linknodesArr(n1-1,n2-1,n3-1,m);}
			if ( n1 ==-1 || n2==-1 || n3 ==-1 ){ ungetc(ch,tre); }
		}
	}
	while (ch!=';');
	return -1;
}
int getcladeArr_UsePartitions(FILE *tre, int whichPartitions, char*** nodeIDsArr_heap){
	int n1, n2, n3;
	char ch;
	n1=-1;
	n2=-1;
	n3=-1;
	do{
		if (specsearchArr_UsePartitions(tre,whichPartitions,nodeIDsArr_heap)==1){
			return tip+(numspecArr[whichPartitions]-1);
		}
		ch = fgetc(tre);
		if (ch==','){
			comma++;
		}
		if (ch==')'){
			if ((ch=fgetc(tre))!=':'){
				ungetc(ch,tre);
			}else{
				do{
					ch=(fgetc(tre));
				}while((ch=='\n')||(ch==' '));
				ungetc(ch,tre);
				fscanf(tre,"%lf",&treeArr[whichPartitions][n3-1].bl);
			}
			return n3;
		}
		if (ch=='('){
			n3=getnodenumbArr(tre);
			n1=getcladeArr_UsePartitions(tre,whichPartitions,nodeIDsArr_heap);
			n2=getcladeArr_UsePartitions(tre,whichPartitions,nodeIDsArr_heap);
			if ( n1!=-1 && n2!=-1 && n3!=-1){linknodesArr(n1-1,n2-1,n3-1,whichPartitions);}
			if ( n1 ==-1 || n2==-1 || n3 ==-1 ){ ungetc(ch,tre); }
		}
	}
	while (ch!=';');
	return -1;
}
void linknodesArr(int i,int j,int node,struct masterArr *m){
	m->tree[node].up[0]=j;
	m->tree[node].up[1]=i;
	m->tree[i].down=node;
	m->tree[j].down=node;
}
int getnodenumbArr(FILE *tre){
	char c;
	int i,j=0;
	fpos_t position;
	i=0;
	fgetpos(tre, &position);
	do{
		c=fgetc(tre);
		if (c==',')
			i++;
		if (c=='(')
			j=j-1;
		if (c==')')
			j++;
	}
	while ((j<0)&&(c!=EOF));
	fsetpos(tre,&position);
	return (i+comma+1);
}
int specsearchArr(FILE *tre, struct masterArr *m, int max_node_name){
	char ch;
	int i=0;
	char specname[max_node_name+1];
	ch = fgetc(tre);
	if ((ch!=')')&&(ch!='(')&&(ch!=',')&&(ch!=' ')&&(ch!='\t')&&(ch!='\n')&&(ch!=EOF)){
		ungetc(ch, tre);
		while ((ch=fgetc(tre))!=':'&&(i<max_node_name)){
			if ( isalpha(ch) || isdigit(ch) || ch=='.' || ch=='_' || ch=='/' || ch=='-' || ch=='|' || ch=='?' || ch=='*' || ch=='&' || ch=='+' || ch=='#' || ch=='\'' || ch=='=' || ch==',' || ch==']' || ch=='[' || ch=='>' || ch=='<'){
				specname[i]=ch;
				i++;
			}
		}
		specname[i]='\0';
		int tmp1=0;
		int found=0;
		for(tmp1=0;tmp1<m->numspec;tmp1++){
			if ( !strcmp(m->names[tmp1],specname) ){
				tip = tmp1+1;
				found=1;
			}
		}
		if (!found) fprintf(stderr, "WARNING: leaf '%s' not found in MSA\n", specname);
		strcpy(m->tree[tip+m->numspec-2].name,specname);
		m->tree[tip+m->numspec-2].up[0]=-1;
		m->tree[tip+m->numspec-2].up[1]=-1;
		fscanf(tre,"%lf",&m->tree[tip+m->numspec-2].bl);
		return 1;
	}
	else {
		ungetc(ch, tre);
		return 0;
	}
}
int specsearchArr_UsePartitions(FILE *tre, int whichPartitions, char*** nodeIDsArr_heap){
	char ch;
	int i=0;
	char specname[30];
	ch = fgetc(tre);
	if ((ch!=')')&&(ch!='(')&&(ch!=',')&&(ch!=' ')&&(ch!='\t')&&(ch!='\n')&&(ch!=EOF)){
		ungetc(ch, tre);
		while ((ch=fgetc(tre))!=':'&&(i<30)){
			if ( isalpha(ch) || isdigit(ch) || ch=='.' || ch=='_' || ch=='/' || ch=='-' || ch=='|' || ch=='?' || ch=='*' || ch=='&' || ch=='+' || ch=='#' || ch=='\'' ){
				specname[i]=ch;
				i++;
			}
		}
		specname[i]='\0';
		int tmp1=0;
		for(tmp1=0;tmp1<numspecArr[whichPartitions];tmp1++){
			if ( !strcmp(nodeIDsArr_heap[whichPartitions][tmp1],specname) ){
				tip = tmp1+1;
			}
		}
		strcpy(treeArr[whichPartitions][tip+numspecArr[whichPartitions]-2].name,specname);
		treeArr[whichPartitions][tip+numspecArr[whichPartitions]-2].up[0]=-1;
		treeArr[whichPartitions][tip+numspecArr[whichPartitions]-2].up[1]=-1;
		fscanf(tre,"%lf",&treeArr[whichPartitions][tip+numspecArr[whichPartitions]-2].bl);
		return 1;
	}
	else {
		ungetc(ch, tre);
		return 0;
	}
}
