/* The SAM parse of tronko-assign: from the aligner's records of one batch, fill each read's
 * candidate slots (concordant_ and discordant_matches_roots/nodes, and with -e the CIGAR and
 * start arrays) in the batch's bwaMatches. Moved unchanged from the vendored BWA's fastmap.c
 * (main_mem's pipeline step 2), including the candidate-slot overflow fix (every slot write
 * checks k < MAX_NUM_BWA_MATCHES) and the read-name buffer of max_readname_length bytes.
 *
 * sam[i] is the SAM text (records ending in newlines, no header) of the i-th read handed to the
 * aligner: for paired input, sam[2k] holds the records of pair k's read 1 and sam[2k+1] those of
 * its read 2. The strings are modified (strtok_r) and not freed. Reads the Tronko fields QNAME,
 * FLAG, RNAME and RNEXT, and POS and CIGAR for -e. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "hashmap.h"
#include "hashmap_base.h"
#include "sam_parse.h"

void sam_parse_batch(const sam_parse_ctx_t *ctx, int n_seqs, char **sam)
{
	int i;
		int j=0;
		int l=0;
		int k=0;
		int no_add=0;
		int success=1;
	HASHMAP(char, leafMap) map;
	hashmap_init(&map, hashmap_hash_string, strcmp);
	for(i=0; i<ctx->ntree; i++){
		for(j=numspecArr[i]-1; j<2*numspecArr[i]-1; j++){
			struct leafMap *l;
			l = malloc(sizeof(*l));
			l->name = treeArr[i][j].name;
			l->root = i;
			l->node = j;
			hashmap_put(&map,l->name,l);
		}
	}
	j=0;
		for (i = 0; i < n_seqs; ++i) {
			/*if (sam[i] && ctx->concordant==1){
				err_fputs(sam[i], stdout);
				char* token;
				char readname[MAXREADNAME];
				char read1[MAX_NODENAME];
				char read2[MAX_NODENAME];
				token=strtok(sam[i],"\n");
				while(success !=0 && token != NULL){
					success = sscanf(token, "%s %*d %s %*d %*d %*s %s %*d %*d %*s %*s %*s %*s %*s %*s",readname,read1,read2);
					if ( j > 0 && strcmp(readname,ctx->results[j-1].readname)==0 ){
						for ( k=0; k<ctx->ntree; k++){
							if ( ctx->results[j-1].matches[k]==-1 ){
								break;
							}
						}
						for (l=0; l<k; l++){
							if (ctx->results[j-1].matches[l]==hashmap_get(&map,read1)){
								no_add=1;
							}
						}
						if ( no_add==0 ){
							ctx->results[j-1].matches[k] = hashmap_get(&map,read1);
							k++;
						}
						no_add=0;
						for(l=0; l<k; l++){
							if ( ctx->results[j-1].matches[l] == hashmap_get(&map,read2) && strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0){
								no_add=1;
							}
						}
						if (no_add==0){
							ctx->results[j-1].matches[k] = hashmap_get(&map,read2);
						}
						no_add=0;
				}else{
					ctx->results[j].readname=token;
					token=strtok(NULL,"\t");
					token=strtok(NULL,"\t");
					ctx->results[j].matches[0] = token;
					j++;
				}*/
			//}else if ( sam[i] && ctx->concordant==0){
			if ( sam[i] ){
				char *token;
				char readname[ctx->max_readname_length];
				char read1[ctx->max_acc_name];
				char read2[ctx->max_acc_name];
				char cigar[MAX_CIGAR];
				int start_position;
				char *rest = sam[i];
				int decimal=0;
				int decimal_general=0;
				while(token=strtok_r(rest,"\n",&rest)){
					success = sscanf(token, "%s %d %s %d %*d %s %s %*d %*d %*s %*s %*s %*s %*s %*s %*[^.,;]",readname,&decimal,read1,&start_position,cigar,read2);
					/*decimal_general = dec2bin(decimal);
					if (ctx->paired==1){
						decimal = dec2bin(decimal);
					}*/
					//if decimal==1, first in pair if decimal==0, second in pair
				//if ( j > 0 && strcmp(readname,ctx->results[j-1].readname)==0 ){
				int whichMat=0;
				if (ctx->paired==1 && j>0){
					if (strcmp(readname,pairedQueryMat->forward_name[ctx->startline+j-1])==0){
						whichMat=1;
					}
				}else if (j>0){
					if (strcmp(readname,singleQueryMat->name[ctx->startline+j-1])==0){
						whichMat=1;
					}
				}
				if ( j > 0 && whichMat==1 ){
					for(k=0; k<MAX_NUM_BWA_MATCHES; k++){
						//if ( ctx->results[j-1].concordant_matches[k] == -1 ){
						//	break;
						//}
						//if ( strlen(ctx->results[j-1].concordant_leaf_matches[k])==0 ){
						//	break;
						//}
						if (ctx->results[j-1].concordant_matches_roots[k] == -1){
							break;
						}
					}
					for(l=0; l<k; l++){
						if ( strcmp(read2,"=")==0 ){
							//if ( ctx->results[j-1].concordant_matches[l]==hashmap_get(&map,read1) ){
							//	no_add=1;
							//}
							//if ( strcmp(ctx->results[j-1].concordant_leaf_matches[l],read1)==0 ){
							//	no_add=1;
							//}
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read1);
							if ( ctx->results[j-1].concordant_matches_roots[l] == leaf_map->root && ctx->results[j-1].concordant_matches_nodes[l] == leaf_map->node){
								no_add=1;
							}
						}
					}
					if (ctx->concordant==1 && no_add==0 && strcmp(read2,"=")!=0){ no_add=1;}
					if (no_add==0 && k < MAX_NUM_BWA_MATCHES && strcmp(read2,"=")==0){
							struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read1);
						ctx->results[j-1].concordant_matches_roots[k] = leaf_map->root;
						ctx->results[j-1].concordant_matches_nodes[k] = leaf_map->node;
						//ctx->results[j-1].concordant_matches[k] = hashmap_get(&map,read1);
						//strcpy(ctx->results[j-1].concordant_leaf_matches[k],read1);
						if ( ctx->results[j-1].use_portion == 1){
							strcpy(ctx->results[j-1].cigars_forward[k],cigar);
							ctx->results[j-1].starts_forward[k] = start_position;
						}
						k++;
					}
					no_add=0;
					for(l=0; l<k; l++){
						//if ( ctx->results[j-1].concordant_matches[l] == hashmap_get(&map,read2) || strcmp(read2,"=") == 0 || strcmp(read2,"*") == 0){
						//	no_add=1;
						//}
						//if ( strcmp(ctx->results[j-1].concordant_leaf_matches[l],read2)==0 || strcmp(read2,"=") == 0 || strcmp(read2,"*") == 0){
						//	no_add=1;
						//}
						if ( strcmp(read2,"=")== 0 || strcmp(read2,"*")==0){
							no_add=1;
							break;
						}
						struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read2);
						if ( ctx->results[j-1].concordant_matches_roots[l] == leaf_map->root && ctx->results[j-1].concordant_matches_nodes[l] == leaf_map->node ){
							no_add=1;
						}
					}
					//if (ctx->concordant==1){ no_add=1;}
					if (k==0){ no_add=1; }
					if (no_add==0 && k < MAX_NUM_BWA_MATCHES && strcmp(read2,"*")!=0 && strcmp(read1,"=")!=0){
						struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read2);
						ctx->results[j-1].concordant_matches_roots[k] = leaf_map->root;
						ctx->results[j-1].concordant_matches_nodes[k] = leaf_map->node;
						//ctx->results[j-1].concordant_matches[k] = hashmap_get(&map,read2);
						//strcpy(ctx->results[j-1].concordant_leaf_matches[k],read2);
						if (ctx->results[j-1].use_portion==1){
							strcpy(ctx->results[j-1].cigars_reverse[k-1],cigar);
							ctx->results[j-1].starts_reverse[k-1] = start_position;
						}
					}
					no_add=0;
					for(k=0; k<MAX_NUM_BWA_MATCHES; k++){
						//if ( ctx->results[j-1].discordant_matches[k]== -1 ){
						//	break;
						//}
						//if ( strlen(ctx->results[j-1].discordant_leaf_matches[k])==0 ){
						//	break;
						//}
						if ( ctx->results[j-1].discordant_matches_roots[k]==-1){
							break;
						}
					}
					for(l=0; l<k; l++){
						if (strcmp(read2,"=") != 0 ){
							//if (ctx->results[j-1].discordant_matches[l]==hashmap_get(&map,read1)){
							//	no_add=1;
							//}
							//if (strcmp(ctx->results[j-1].discordant_leaf_matches[l],read1)==0){
							//	no_add=1;
							//}
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read1);
							if ( ctx->results[j-1].discordant_matches_roots[l] == leaf_map->root && ctx->results[j-1].discordant_matches_nodes[l] == leaf_map->node ){
								no_add=1;
							}
						}
					}
					if (no_add==0 && k < MAX_NUM_BWA_MATCHES && strcmp(read2,"=") != 0){
							struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read1);
						ctx->results[j-1].discordant_matches_roots[k] = leaf_map->root;
						ctx->results[j-1].discordant_matches_nodes[k] = leaf_map->node;
						//ctx->results[j-1].discordant_matches[k]=hashmap_get(&map,read1);
						//strcpy(ctx->results[j-1].discordant_leaf_matches[k],read1);
						if (ctx->results[j-1].use_portion==1 && decimal == 1){
							strcpy(ctx->results[j-1].cigars_forward[k],cigar);
							ctx->results[j-1].starts_forward[k] = start_position;
						}
						k++;
					}
					no_add=0;
					for(l=0; l<k; l++){
						//if (ctx->results[j-1].discordant_matches[l] == hashmap_get(&map,read2) || strcmp(read2,"=") ==0 || strcmp(read2,"*") == 0){
						//	no_add=1;
						//}
						if (strcmp(read2,"*")==0 || strcmp(read2,"=")==0){
							no_add=1;
							break;
						}
						struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read2);
						if ( ctx->results[j-1].discordant_matches_roots[l] == leaf_map->root && ctx->results[j-1].discordant_matches_nodes[l] == leaf_map->node){
							no_add=1;
						}
						//if ( strcmp(ctx->results[j-1].discordant_leaf_matches[l],read2)==0 || strcmp(read2,"=")==0 || strcmp(read2,"*") == 0 ){
						//	no_add=1;
						//}
					}
					if (no_add==0 && k < MAX_NUM_BWA_MATCHES && strcmp(read2,"=")!=0 && strcmp(read2,"*")!=0){
						struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read2);
						ctx->results[j-1].discordant_matches_roots[k] = leaf_map->root;
						ctx->results[j-1].discordant_matches_nodes[k] = leaf_map->node;
						//ctx->results[j-1].discordant_matches[k] = hashmap_get(&map,read2);
						//strcpy(ctx->results[j-1].discordant_leaf_matches[k],read2);
						if (ctx->results[j-1].use_portion==1){
							strcpy(ctx->results[j-1].cigars_reverse[k],cigar);
							ctx->results[j-1].starts_reverse[k] = start_position;
						}
					}
					no_add=0;
					if ( decimal == 0 && k < MAX_NUM_BWA_MATCHES && ctx->results[j-1].use_portion==1){
						strcpy(ctx->results[j-1].cigars_reverse[k],cigar);
						ctx->results[j-1].starts_reverse[k] = start_position;
					}
				}else if (j==0 && strcmp(read2,"=")==0){
					//strcpy(ctx->results[j].readname,readname);
							struct leafMap *leaf_map;
					leaf_map=hashmap_get(&map,read1);
					ctx->results[j].concordant_matches_roots[0] = leaf_map->root;
					ctx->results[j].concordant_matches_nodes[0] = leaf_map->node;
					//ctx->results[j].concordant_matches[0] = hashmap_get(&map,read1);
					//strcpy(ctx->results[j].concordant_leaf_matches[0],read1);
					if (ctx->results[j].use_portion==1){
						strcpy(ctx->results[j].cigars_forward[0],cigar);
						ctx->results[j].starts_forward[0] = start_position;
					}
					//if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && ctx->results[j].concordant_matches[0] != hashmap_get(&map,read2)){
					//if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && strcmp(ctx->results[j].concordant_leaf_matches[0],read2) != 0){
					if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && (ctx->results[j].concordant_matches_roots[0] != leaf_map->root && ctx->results[j].concordant_matches_nodes[0] != leaf_map->node)){
							struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read2);
						ctx->results[j].concordant_matches_roots[1] = leaf_map->root;
						ctx->results[j].concordant_matches_nodes[1] = leaf_map->node;
						//ctx->results[j].concordant_matches[1] = hashmap_get(&map,read2);
						//strcpy(ctx->results[j].concordant_leaf_matches[1],read2);
						if (ctx->results[j].use_portion==1){
							strcpy(ctx->results[j].cigars_reverse[0],cigar);
							ctx->results[j].starts_reverse[0] = start_position;
						}
					}
					j++;
				}else if (j==0 && strcmp(read2,"=")!=0 && ctx->paired != 0 && strcmp(read1,"*")!=0){
					//strcpy(ctx->results[j].readname,readname);
					struct leafMap *leaf_map;
					leaf_map=hashmap_get(&map,read1);
					ctx->results[j].discordant_matches_roots[0] = leaf_map->root;
					ctx->results[j].discordant_matches_nodes[0] = leaf_map->node;
					//ctx->results[j].discordant_matches[0] = hashmap_get(&map,read1);
					//strcpy(ctx->results[j].discordant_leaf_matches[0],read1); 
					if (ctx->results[j].use_portion==1){
						strcpy(ctx->results[j].cigars_forward[0],cigar);
						ctx->results[j].starts_forward[0] = start_position;
					}
					//if (strcmp(read2,"*")!=0 && ctx->results[j].discordant_matches[0] != hashmap_get(&map,read2)){
					//if ( strcmp(read2,"*") != 0 && strcmp(ctx->results[j].discordant_leaf_matches[0],read2) != 0 ){
					if ( strcmp(read2,"*") != 0 && (ctx->results[j].discordant_matches_roots[0] != leaf_map->root && ctx->results[j].discordant_matches_nodes[0] != leaf_map->node)){
							struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read2);
						ctx->results[j].discordant_matches_roots[1] = leaf_map->root;
						ctx->results[j].discordant_matches_nodes[1] = leaf_map->node;
						//ctx->results[j].discordant_matches[1] = hashmap_get(&map,read2);
						//strcpy(ctx->results[j].discordant_leaf_matches[1],read2);
						if (ctx->results[j].use_portion==1){
							strcpy(ctx->results[j].cigars_reverse[0],cigar);
							ctx->results[j].starts_reverse[0] = start_position;
						}
					}
					j++;
				}else{
					if (ctx->paired != 0 && strcmp(pairedQueryMat->forward_name[ctx->startline+j],readname)==0 && strcmp(read2,"=")==0){
						if ( decimal == 1 || decimal==2){
							//strcpy(ctx->results[j].readname,readname);
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read1);
							ctx->results[j].concordant_matches_roots[0] = leaf_map->root;
							ctx->results[j].concordant_matches_nodes[0] = leaf_map->node;
							//ctx->results[j].concordant_matches[0] = hashmap_get(&map,read1);
							//strcpy(ctx->results[j].concordant_leaf_matches[0],read1);
							if (ctx->results[j].use_portion==1){
								strcpy(ctx->results[j].cigars_forward[0],cigar);
								ctx->results[j].starts_forward[0] = start_position;
							}
						}
						//if (ctx->results[j].use_portion==1){
						//	strcpy(ctx->results[j].cigars_reverse[0],cigar);
						//	ctx->results[j].starts_reverse[0] = start_position;
						//}
						//if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && ctx->results[j].concordant_matches[0] != hashmap_get(&map,read2)){
						//if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && strcmp(ctx->results[j].concordant_leaf_matches[0],read2) != 0 ){
						if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && strcmp(read1,read2) != 0 ){
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read2);
							ctx->results[j].concordant_matches_roots[1] = leaf_map->root;
							ctx->results[j].concordant_matches_nodes[1] = leaf_map->node;
							//ctx->results[j].concordant_matches[1] = hashmap_get(&map,read2);
							//strcpy(ctx->results[j].concordant_leaf_matches[1],read2);
							if (ctx->results[j].use_portion==1){
								strcpy(ctx->results[j].cigars_reverse[0],cigar);
								ctx->results[j].starts_reverse[0] = start_position;
							}
						}
						j++;
					}else if (ctx-> paired != 0 && strcmp(pairedQueryMat->forward_name[ctx->startline+j],readname)==0 && strcmp(read2,"=")!=0 && strcmp(read1,"*")!=0){
						//strcpy(ctx->results[j].readname,readname);
							struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read1);
						ctx->results[j].discordant_matches_roots[0] = leaf_map->root;
						ctx->results[j].discordant_matches_nodes[0] = leaf_map->node;
						//ctx->results[j].discordant_matches[0] = hashmap_get(&map,read1);
						//if (strcmp(read2,"*") != 0 && ctx->results[j].discordant_matches[0] != hashmap_get(&map,read2) ){
						//if (strcmp(read2,"*") != 0 && strcmp(ctx->results[j].discordant_leaf_matches[0],read2) != 0 ){
						if (strcmp(read2,"*") != 0 && strcmp(read1,read2)!=0 ){
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read2);
							ctx->results[j].discordant_matches_roots[1] = leaf_map->root;
							ctx->results[j].discordant_matches_nodes[1] = leaf_map->node;
							//ctx->results[j].discordant_matches[1] = hashmap_get(&map,read2);
							//strcpy(ctx->results[j].discordant_leaf_matches[1],read2);
							if (ctx->results[j].use_portion==1){
								strcpy(ctx->results[j].cigars_reverse[0],cigar);
								ctx->results[j].starts_reverse[0] = start_position;
							}
						}
						j++;
					}else if (ctx->paired==0 && strcmp(singleQueryMat->name[ctx->startline+j],readname)==0 && strcmp(read2,"=")==0){
						//strcpy(ctx->results[j].readname,readname);
							struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read1);
						ctx->results[j].concordant_matches_roots[0] = leaf_map->root;
						ctx->results[j].concordant_matches_nodes[0] = leaf_map->node;
						//ctx->results[j].concordant_matches[0] = hashmap_get(&map,read1);
						//strcpy(ctx->results[j].concordant_leaf_matches[0],read1);
						if (ctx->results[j].use_portion==1){
							strcpy(ctx->results[j].cigars_forward[0],cigar);
							ctx->results[j].starts_forward[0] = start_position;
						}
						//if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && strcmp(ctx->results[j].concordant_leaf_matches[0],read2) != 0 ){
						if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && strcmp(read1,read2)!=0 ){
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read2);
							ctx->results[j].concordant_matches_roots[1] = leaf_map->root;
							ctx->results[j].concordant_matches_nodes[1] = leaf_map->node;
							//ctx->results[j].concordant_matches[1] = hashmap_get(&map,read2);
							//strcpy(ctx->results[j].concordant_leaf_matches[1],read2);
							if (ctx->results[j].use_portion==1){
								strcpy(ctx->results[j].cigars_reverse[0],cigar);
								ctx->results[j].starts_reverse[0] = start_position;
							}
						}
						j++;
					}else if (ctx->paired==0 && strcmp(singleQueryMat->name[ctx->startline+j],readname)==0 && strcmp(read2,"=")!=0){
						//strcpy(ctx->results[j].readname,readname);
						if ( strcmp(read1,"*") != 0 ){
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read1);
							ctx->results[j].discordant_matches_roots[0] = leaf_map->root;
							ctx->results[j].discordant_matches_nodes[0] = leaf_map->node;
							//ctx->results[j].discordant_matches[0] = hashmap_get(&map,read1);
							//strcpy(ctx->results[j].discordant_leaf_matches[0],read1);
							if(ctx->results[j].use_portion==1){
								strcpy(ctx->results[j].cigars_forward[0],cigar);
								ctx->results[j].starts_forward[0]=start_position;
							}
							if (strcmp(read2,"*") != 0 && strcmp(read1,read2) != 0 ){
							struct leafMap *leaf_map;
								leaf_map=hashmap_get(&map,read2);
								ctx->results[j].discordant_matches_roots[1] = leaf_map->root;
								ctx->results[j].discordant_matches_nodes[1] = leaf_map->node;
								//ctx->results[j].discordant_matches[1] = hashmap_get(&map,read2);
								//strcpy(ctx->results[j].discordant_leaf_matches[1],read2);
								if (ctx->results[j].use_portion==1){
									strcpy(ctx->results[j].cigars_reverse[0],cigar);
									ctx->results[j].starts_reverse[0] = start_position;
								}
							}
						}else{
							ctx->results[j].concordant_matches_roots[0]=-1;
							ctx->results[j].concordant_matches_nodes[0]=-1;
							ctx->results[j].discordant_matches_nodes[0]=-1;
							ctx->results[j].discordant_matches_roots[0]=-1;
							//ctx->results[j].concordant_matches[0]=-2;
							//ctx->results[j].discordant_matches[0]=-2;
						}
						j++;	
					}else{
						if (ctx->paired != 0){
							while(strcmp(pairedQueryMat->forward_name[ctx->startline+j],readname)!=0){
								//strcpy(ctx->results[j].readname,pairedQueryMat->forward_name[ctx->startline+j]);
								ctx->results[j].concordant_matches_roots[0]=-1;
								ctx->results[j].concordant_matches_nodes[0]=-1;
								ctx->results[j].discordant_matches_nodes[0]=-1;
								ctx->results[j].discordant_matches_roots[0]=-1;
								//ctx->results[j].concordant_matches[0]=-2;
								//ctx->results[j].discordant_matches[0]=-2;
								j++;
							}
						}else{
							while(strcmp(singleQueryMat->name[ctx->startline+j],readname)!=0){
								//strcpy(ctx->results[j].readname,singleQueryMat->name[ctx->startline+j]);
								ctx->results[j].concordant_matches_roots[0]=-1;
								ctx->results[j].concordant_matches_nodes[0]=-1;
								ctx->results[j].discordant_matches_nodes[0]=-1;
								ctx->results[j].discordant_matches_roots[0]=-1;
								//ctx->results[j].concordant_matches[0]=-2;
								//ctx->results[j].discordant_matches[0]=-2;
								j++;
							}
						}
					}
				}
				}
			}
		}
	struct leafMap *blob;
	hashmap_foreach_data(blob,&map){
		free(blob);
	}
	hashmap_cleanup(&map);
}
