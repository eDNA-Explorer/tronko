#include "placement.h"
#include <math.h>

int perform_WFA_alignment(cigar_t* const cigar, mm_allocator_t* mm_allocator,char* seq1, char* seq2,char* const pattern_alg,char* const text_alg, char* const ops_alg, int begin_offset, int end_offset){
	char* const operations = cigar->operations;
	int k, alg_pos =0, pattern_pos= 0, text_pos =0;
	for (k=begin_offset;k<end_offset;++k) {
		switch (operations[k]) {
			case 'M':
				if (seq1[pattern_pos] != seq2[text_pos]) {
					pattern_alg[alg_pos] = seq1[pattern_pos];
					ops_alg[alg_pos] = 'X';
					text_alg[alg_pos++] = seq2[text_pos];
				}else{
					pattern_alg[alg_pos] = seq1[pattern_pos];
					ops_alg[alg_pos] = '|';
					text_alg[alg_pos++] = seq2[text_pos];
				}
				pattern_pos++; text_pos++;
				break;
			case 'X':
				if (seq1[pattern_pos] != seq2[text_pos]) {
					pattern_alg[alg_pos] = seq1[pattern_pos++];
					ops_alg[alg_pos] = ' ';
					text_alg[alg_pos++] = seq2[text_pos++];
				}else{
					pattern_alg[alg_pos] = seq1[pattern_pos++];
					ops_alg[alg_pos] = 'X';
					text_alg[alg_pos++] = seq2[text_pos++];
				}
				break;
			case 'I':
				pattern_alg[alg_pos] = '-';
				ops_alg[alg_pos] = ' ';
				text_alg[alg_pos++] = seq2[text_pos++];
				break;
			case 'D':
				pattern_alg[alg_pos] = seq1[pattern_pos++];
				ops_alg[alg_pos] = ' ';
				text_alg[alg_pos++] = '-';
				break;
			default:
				break;
		}
	}
	k=0;
	while (pattern_pos < strlen(seq1)) {
		pattern_alg[alg_pos+k] = seq1[pattern_pos++];
		ops_alg[alg_pos+k] = '?';
		++k;
	}
	while (text_pos < strlen(seq2)) {
		text_alg[alg_pos+k] = seq2[text_pos++];
		ops_alg[alg_pos+k] = '?';
		++k;
	}
	int alignment_length = strlen(pattern_alg);
	return alignment_length;
}

void place_paired( char *query_1, char *query_2, char **rootSeqs, int numberOfTotalRoots, int *positions, char *locQuery, type_of_PP ***nodeScores, double **voteRoot, int number_of_matches , int **leaf_coordinates, int paired, type_of_PP* minimum_score, char *alignments_dir, char *forward_name, char *reverse_name, int print_alignments, char *leaf_sequence, int *positionsInRoot, int maxNumSpec, int* starts_forward, char** cigars_forward, int* starts_reverse, char** cigars_reverse, int print_alignments_to_file, int use_leaf_portion, int padding, int max_query_length, int max_numbase, int print_all_nodes, int early_termination, type_of_PP strike_box, int max_strikes, int enable_pruning, type_of_PP pruning_threshold, type_of_PP best_leaf_threshold, int best_leaf_max_votes){
	int i, j, k, node, match;
	type_of_PP forward_mismatch, reverse_mismatch;
	forward_mismatch=0;
	reverse_mismatch=0;
	/*affine_penalties.match = 0;
	affine_penalties.mismatch =4;
	affine_penalties.gap_opening=6;
	affine_penalties.gap_extension = 2;
	#define BUFFER_SIZE_8M   (1ul<<23)*/
	wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
	attributes.distance_metric = gap_affine;
	attributes.affine_penalties.mismatch =4;
	attributes.affine_penalties.gap_opening = 6;
	attributes.affine_penalties.gap_extension = 2;
	attributes.alignment_form.span = alignment_endsfree;
	/* Create WFA aligner once and reuse across all matches */
	wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
	for(match=0; match<number_of_matches; match++){
		{
			int clear_len = use_leaf_portion ? (max_query_length+max_query_length+2*padding+1) : (max_query_length+max_numbase+1);
			memset(positions, -1, clear_len * sizeof(int));
			memset(locQuery, '\0', clear_len);
			memset(leaf_sequence, '\0', clear_len);
			memset(positionsInRoot, -1, clear_len * sizeof(int));
		}
		int query_length = strlen(query_1);
		if (leaf_coordinates[match][0] == -1 || leaf_coordinates[match][1] == -1) {
			continue;  /* Skip invalid leaf coordinates */
		}
		if (use_leaf_portion==1){
			if ( cigars_forward[match][0] == '*'){ break; }
			if ( cigars_forward[match][0] == '\0'){ break; }
		}
		if (use_leaf_portion == 1 && starts_forward[match] != -1){
			int start_position = getStartPosition(starts_forward[match],leaf_coordinates[match][0],leaf_coordinates[match][1],padding);
			int end_position = getEndPosition(cigars_forward[match],leaf_coordinates[match][0],leaf_coordinates[match][1],start_position+padding,padding);
			getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,start_position,end_position);
		}else{
			getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,0,numbaseArr[leaf_coordinates[match][0]]);
		}
		int leaf_length = strlen(leaf_sequence);
		if (leaf_length > 0){
		/*mm_allocator_t* mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
		affine_wavefronts_t* affine_wavefronts;
		char* const pattern_alg=mm_allocator_calloc(mm_allocator,leaf_length+query_length+1,char,true);
		char* const text_alg=mm_allocator_calloc(mm_allocator,leaf_length+query_length+1,char,true);*/
		// Align using reused WFA aligner
		wavefront_align(wf_aligner,leaf_sequence,leaf_length,query_1,query_length);
		char* const pattern_alg = mm_allocator_calloc(wf_aligner->mm_allocator,leaf_length+query_length+1,char,true);
		char* const ops_alg = mm_allocator_calloc(wf_aligner->mm_allocator,leaf_length+query_length+1,char,true);
		char* const text_alg = mm_allocator_calloc(wf_aligner->mm_allocator,leaf_length+query_length+1,char,true);
		int alignment_length = 0;
		/*affine_wavefronts = affine_wavefronts_new_complete(leaf_length,query_length,&affine_penalties,NULL,mm_allocator);
		affine_wavefronts_align(affine_wavefronts,leaf_sequence,leaf_length,query_1,query_length);*/
		//pattern_alg = mm_allocator_calloc(mm_allocator,leaf_length+query_length+1,char,true);
		//text_alg = mm_allocator_calloc(mm_allocator,leaf_length+query_length+1,char,true);
		alignment_length = perform_WFA_alignment(wf_aligner->cigar,wf_aligner->mm_allocator,leaf_sequence,query_1,pattern_alg,text_alg,ops_alg,wf_aligner->cigar->begin_offset,wf_aligner->cigar->end_offset);
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		int alength=0;
		int keepTrackOfPosInRoot=0;
		int loop=0;
		int lengthOfResultA=0;
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (aln->result_a[i]=='-'){*/
				/*positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;*/
				//keepTrackOfPosInRoot++;
			/*}else if (aln->result_a[i]=='N'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;
				keepTrackOfPosInRoot++;
			}else{
				if (alength != 0 && aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
				//if ( aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
					locQuery[alength]=aln->result_b[loop];
					positions[alength]=keepTrackOfPosInRoot;
					alength++;
				}
				if ( aln->result_b[loop] == '-' && aln->result_a[loop] != '-'){
					keepTrackOfPosInRoot++;
					loop++;
				}else{
					if (aln->result_a[i]!='-'){
						locQuery[alength]=aln->result_b[loop];
						positions[alength]=keepTrackOfPosInRoot;
						keepTrackOfPosInRoot++;
						alength++;
						loop++;
						lengthOfResultA++;
					}
				}
			}
		}*/
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (keepTrackOfPosInRoot > leaf_length){
				positions[0]=-1;
				break;
			}
			//if (aln->result_a[i] != '-' && aln->result_b[i] != '-'){
			if (aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				if (aln->result_a[i] != aln->result_b[i] && match==0){
					forward_mismatch++;
				}
			}
			if (alength != 0 && lengthOfResultA<query_length && aln->result_b[i] =='-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				keepTrackOfPosInRoot++;
				if (aln->result_a[i] != aln->result_b[i] && match==0 && aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
					forward_mismatch++;
				}
			}else if (aln->result_a[i] != '-'){
			//}else if (aln->result_a[i] != 'N'){
				keepTrackOfPosInRoot++;
			}
		}
		locQuery[alength]='\0';*/
		int counterInRoot=0;
		int counterinB=0;
			for (i=0; pattern_alg[i] != '\0'; i++){
				if (pattern_alg[i] != '-' && text_alg[i] != '-' ){
					positions[alength] = positionsInRoot[counterInRoot];
					locQuery[alength] = text_alg[i];
					alength++;
					if (pattern_alg[i] != text_alg[i] && match==0){
						forward_mismatch++;
					}
				}
				if ( pattern_alg[i] == '-' && text_alg[i] != '-' && alength > 0 && counterinB<query_length){
					positions[alength]=-1;
					locQuery[alength] = text_alg[i];
					alength++;
				}
				if (pattern_alg[i] != '-' && text_alg[i] == '-' && alength > 0 && counterinB<query_length){
					positions[alength] = positionsInRoot[counterInRoot];
					locQuery[alength] = text_alg[i];
					alength++;
				}
				if ( pattern_alg[i] != '-'){
					counterInRoot++;
				}
				if ( text_alg[i] != '-'){
					counterinB++;
				}
			}
		locQuery[alength]='\0';
		//if (access(alignmentFileName, F_OK ) != -1 && match==0){
		//	printToFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, forward_name);
		//}else if (match==0){
		//	createNewFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, forward_name);
		//}
		if (match==0 && print_alignments_to_file==1){
			char alignmentFileName[BUFFER_SIZE];
			snprintf(alignmentFileName,BUFFER_SIZE,"%s/%s.fasta",alignments_dir,treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			snprintf(alignmentFileName,BUFFER_SIZE,"%s",alignments_dir);
			printToFile_WFA(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, pattern_alg, text_alg, forward_name, query_length, positionsInRoot);
		}
		if (print_alignments==1){
			int breaks = 100;
			i=0;
			while(pattern_alg[i] != '\0'){
			printf("\n");
			printf("%s\n",forward_name);
			  printf("query_1\t\t\t\t\t");
			  for(i=breaks-100; i<breaks; i++){
				  if (pattern_alg[i]=='\0'){ break; }
					printf("%c",text_alg[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=breaks-100; i<breaks; i++){
				  if (pattern_alg[i]=='\0'){ break;}
					  printf("%c",pattern_alg[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=breaks-100; i<breaks; i++){
				  if (pattern_alg[i]=='\0'){break;}
					if ( pattern_alg[i] == text_alg[i] && pattern_alg[i] != '-' && text_alg[i] != '-' ){
						printf("|");
					}else if ( pattern_alg[i] != text_alg[i] && pattern_alg[i] != 'N' && text_alg[i] != '-' && pattern_alg[i] != '-'){
						printf("*");
					}else{
						printf(" ");
					}
			  }
			  printf("\n\n");
			  breaks += 100;
			}
			  /*printf("query_1\t\t\t\t\t");
			  for(i=101;i<200;i++){
					printf("%c",text_alg[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=101;i<200;i++){
					printf("%c",pattern_alg[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=101;i<200;i++){
					  if ( pattern_alg[i] == text_alg[i] && pattern_alg[i] != '-' && text_alg[i] != '-' ){
						  printf("|");
					  }else if ( pattern_alg[i] != text_alg[i] && pattern_alg[i] != 'N' && text_alg[i] != '-'){
						  printf("*");
					  }else{
						  printf(" ");
					  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=201; i<300; i++){
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=201;i<300;i++){
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=201;i<300;i++){
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
			  printf("*");
			  }else{
			  printf(" ");
			  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=301; i<400 ; i++){
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=301;i<400; i++){
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=301;i<400;i++){
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
				printf("*");
			  }else{
				  printf(" ");
			  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=401; i<500 ; i++){
				printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=401; i<500 ; i++){
				  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=401; i<500 ; i++){
				  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
						  }else{
							  printf(" ");
						  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=501; i<600; i++){
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=501;i<600; i++){
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=501;i<600;i++){
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
				  printf("*");
			  }else{
				  printf(" ");
			  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=601; i<700; i++){
				  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=601;i<700; i++){
				  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=601;i<700;i++){
				  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
						  }else{
							  printf(" ");
						  }
			  }
			  printf("\n\n");
			printf("query_1\t\t\t\t\t");
			for(i=701; i<800; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=701;i<800; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=701;i<800;i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
					}else{
						printf(" ");
					}
			}
			printf("\n\n");
			printf("query_1\t\t\t\t\t");
			for(i=801;i<900; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=801;i<900; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=801;i<900;i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					 printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						   printf("*");
					  }else{
						  printf(" ");
					  }
			}
			 printf("\n\n");
			printf("query_1\t\t\t\t\t");
			for(i=901; aln->result_b[i] != '\0'; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=901;aln->result_b[i]!='\0'; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=901;aln->result_b[i]!='\0';i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
					}else{
						printf(" ");
					}
			}
			printf("\n\n");*/
		}
		//assignScores_Arr_paired(topRoots[rootNum],rootArr[topRoots[rootNum]],locQuery, positions, nodeScores, alength);
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		//printf("assignScores_Arr_paired...\n");
		FILE* site_scores_file;
		if ( print_all_nodes == 1 && match > 0 && access("site_scores.txt", F_OK ) != -1 ){
			if (( site_scores_file = fopen("site_scores.txt","a")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
			//fprintf(site_scores_file,"%s\t",forward_name);
		}
		if ( print_all_nodes == 1 && match > 0 ){
			if (( site_scores_file = fopen("site_scores.txt","w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
			fprintf(site_scores_file,"Readname\tTree_Number\tNode_Number\tSite\tScore\n");
			//fprintf(site_scores_file,"%s\t",forward_name);
		}
		// Initialize early termination state for this tree
		type_of_PP best_score = -9999999999999999;
		int strikes = 0;
		type_of_PP strike_box_threshold = Cinterval * strike_box;
		type_of_PP pruning_threshold_calc = Cinterval * pruning_threshold;

		assignScores_Arr_paired(leaf_coordinates[match][0],rootArr[leaf_coordinates[match][0]],locQuery, positions, nodeScores, alength, match,print_all_nodes,site_scores_file,forward_name,
		    early_termination, &best_score, &strikes, strike_box_threshold, max_strikes,
		    enable_pruning, pruning_threshold_calc);
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		if ( print_all_nodes == 1){
			fclose(site_scores_file);
		}
		//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		mm_allocator_free(wf_aligner->mm_allocator,pattern_alg);
		mm_allocator_free(wf_aligner->mm_allocator,ops_alg);
		mm_allocator_free(wf_aligner->mm_allocator,text_alg);
	}
	}
	if (paired==1){
	for(match=0; match<number_of_matches; match++){
		{
			int clear_len = use_leaf_portion ? (max_query_length+max_query_length+2*padding+1) : (max_query_length+max_numbase+1);
			memset(positions, -1, clear_len * sizeof(int));
			memset(locQuery, '\0', clear_len);
			memset(leaf_sequence, '\0', clear_len);
			memset(positionsInRoot, -1, clear_len * sizeof(int));
		}
		int query_length = strlen(query_2);
		if (use_leaf_portion==1){
			if ( cigars_reverse[match][0] == '*'){
				break;
			}
			if ( cigars_reverse[match][0] == '\0'){ break; }
		}
		if (use_leaf_portion==1){
			if (starts_reverse[match]!=1){
		//needleman_wunsch_align(rootSeqs[topRoots[rootNum]], query_2, scoring, nw, aln);
		//char *leaf_sequence = (char *)malloc(numbaseArr[leaf_coordinates[rootNum][0]]*sizeof(char));
		//getSequenceInNode(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence);
				int start_position = getStartPosition(starts_reverse[match],leaf_coordinates[match][0],leaf_coordinates[match][1],padding);
				int end_position = getEndPosition(cigars_reverse[match],leaf_coordinates[match][0],leaf_coordinates[match][1],start_position+padding,padding);
				getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,start_position,end_position);
			}
		}else{
			getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,0,numbaseArr[leaf_coordinates[match][0]]);
		}
		int leaf_length = strlen(leaf_sequence);
		if (leaf_length > 0){
		//printf("leaf sequence: ");
		//for(i=0; i<strlen(leaf_sequence); i++){
		//	printf("%c",leaf_sequence[i]);
		//}
		//printf("\n");
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		//printf("aligning 2nd pair...\n");
		/*mm_allocator_t* mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
		affine_wavefronts_t* affine_wavefronts;
		char* const pattern_alg=mm_allocator_calloc(mm_allocator,leaf_length+query_length+1,char,true);
		char* const text_alg=mm_allocator_calloc(mm_allocator,leaf_length+query_length+1,char,true);*/
		wavefront_align(wf_aligner,leaf_sequence,leaf_length,query_2,query_length);
		char* const pattern_alg = mm_allocator_calloc(wf_aligner->mm_allocator,leaf_length+query_length+1,char,true);
		char* const ops_alg = mm_allocator_calloc(wf_aligner->mm_allocator,leaf_length+query_length+1,char,true);
		char* const text_alg = mm_allocator_calloc(wf_aligner->mm_allocator,leaf_length+query_length+1,char,true);
		int alignment_length=0;
		alignment_length = perform_WFA_alignment(wf_aligner->cigar,wf_aligner->mm_allocator,leaf_sequence,query_2,pattern_alg,text_alg,ops_alg,wf_aligner->cigar->begin_offset,wf_aligner->cigar->end_offset);
			/*affine_wavefronts = affine_wavefronts_new_complete(leaf_length,query_length,&affine_penalties,NULL,mm_allocator);
			affine_wavefronts_align(affine_wavefronts,leaf_sequence,leaf_length,query_2,query_length);
			alignment_length = perform_WFA_alignment(affine_wavefronts,mm_allocator,leaf_sequence,query_2,pattern_alg,text_alg);*/
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		/*printf("leaf sequence: ");
		for(i=0; i<strlen(leaf_sequence); i++){
			printf("%c",leaf_sequence[i]);
		}*/
		int alength=0;
		int keepTrackOfPosInRoot=0;
		int loop=0;
		int lengthOfResultA=0;
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (aln->result_a[i]=='-'){*/
				/*positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;*/
				//keepTrackOfPosInRoot++;
			/*}else if (aln->result_a[i]=='N'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;
				keepTrackOfPosInRoot++;
			}else{
				if (alength != 0 && aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
				//if ( aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
					locQuery[alength]=aln->result_b[loop];
					positions[alength]=keepTrackOfPosInRoot;
					alength++;
				}
				if ( aln->result_b[loop] == '-' && aln->result_a[loop] != '-'){
					keepTrackOfPosInRoot++;
					loop++;
				}else{
					if (aln->result_a[i]!='-'){
						locQuery[alength]=aln->result_b[loop];
						positions[alength]=keepTrackOfPosInRoot;
						keepTrackOfPosInRoot++;
						alength++;
						loop++;
						lengthOfResultA++;
					}
				}
			}
		}*/
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (keepTrackOfPosInRoot > leaf_length){
				positions[0]=-1;
				break;
			}
			//if (aln->result_a[i] != '-' && aln->result_b[i] != '-'){
			if (aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				if (aln->result_a[i] != aln->result_b[i] && match==0){
					reverse_mismatch++;
				}
			}
			if (alength != 0 && lengthOfResultA<query_length && aln->result_b[i]=='-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				keepTrackOfPosInRoot++;
				if (aln->result_a[i] != aln->result_b[i] && match==0 && aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
					reverse_mismatch++;
				}
			}else if (aln->result_a[i] != '-'){
			//}else if (aln->result_a[i] != 'N'){
				keepTrackOfPosInRoot++;
			}
		}
		locQuery[alength]='\0';*/
		int counterInRoot=0;
		int counterinB=0;
			for (i=0; pattern_alg[i] != '\0'; i++){
				if (pattern_alg[i] != '-' && text_alg[i] != '-' ){
					positions[alength] = positionsInRoot[counterInRoot];
					locQuery[alength] = text_alg[i];
					alength++;
					if (pattern_alg[i] != text_alg[i] && match==0){
						reverse_mismatch++;
					}
				}
				if ( pattern_alg[i] == '-' && text_alg[i] != '-' && alength > 0 && counterinB<query_length){
					positions[alength]=-1;
					locQuery[alength] = text_alg[i];
					alength++;
				}
				if (pattern_alg[i] != '-' && text_alg[i] == '-' && alength > 0 && counterinB<query_length){
					positions[alength] = positionsInRoot[counterInRoot];
					locQuery[alength] = text_alg[i];
					alength++;
				}
				if ( pattern_alg[i] != '-'){
					counterInRoot++;
				}
				if ( text_alg[i] != '-'){
					counterinB++;
				}
			}
		locQuery[alength]='\0';
		//if (access(alignmentFileName, F_OK ) != -1 && match==0){
		//	printToFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, reverse);
		//}else if (match==0){
		//	createNewFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, readname);
		//}
		if (match==0 && print_alignments_to_file==1){
			char alignmentFileName[1000];
			snprintf(alignmentFileName,1000,"%s/%s.fasta",alignments_dir,treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			printToFile_WFA(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, pattern_alg, text_alg, reverse_name, query_length, positionsInRoot);
		}
		if (print_alignments==1){
			int breaks = 100;
			i=0;
			while(pattern_alg[i] != '\0'){
			printf("\n");
			printf("%s\n",reverse_name);
			printf("query_2\t\t\t\t\t");
			for(i=breaks-100; i<breaks; i++){
				if (text_alg[i]=='\0'){ break; }
					printf("%c",text_alg[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=breaks-100; i<breaks; i++){
				if (pattern_alg[i]=='\0'){ break; }
					printf("%c",pattern_alg[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=breaks-100; i<breaks; i++){
				if (pattern_alg[i]=='\0'){ break; }
					if ( pattern_alg[i] == text_alg[i] && pattern_alg[i] != '-' && text_alg[i] != '-' ){
						printf("|");
					}else if ( pattern_alg[i] != text_alg[i] && pattern_alg[i] != 'N' && text_alg[i] != '-' && pattern_alg[i] != '-'){
						printf("*");
					}else{
						printf(" ");
					}
			}
			printf("\n\n");
			breaks += 100;
			}
			/*
			printf("query_2\t\t\t\t\t");
			for(i=101;i<200;i++){
			  printf("%c",text_alg[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=101;i<200;i++){
			  printf("%c",pattern_alg[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=101;i<200;i++){
			  if(pattern_alg[i] == text_alg[i] && pattern_alg[i] !='-' && text_alg[i] != '-'){
			  printf("|");
			  }else if ( pattern_alg[i] != text_alg[i] && pattern_alg[i] != 'N' && text_alg[i] != '-' ){
			  printf("*");
			  }else{
			  printf(" ");
			  }
			}
			printf("\n\n");
			  printf("query_2\t\t\t\t\t");
			  for(i=201; i<300; i++){
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=201;i<300;i++){
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=201;i<300;i++){
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
			  printf("*");
			  }else{
			  printf(" ");
			  }
			  }
			  printf("\n\n");
			  printf("\nquery_2\t\t\t\t\t");
			  for(i=301; i<400; i++){
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=301;i<400; i++){
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=301;i<400;i++){
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
				  printf("*");
			  }else{
				  printf(" ");
			  }
			  }
			  printf("\n\n");
			  printf("\nquery_2\t\t\t\t\t");
			  for(i=401; i<500; i++){
				  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=401; i<500; i++){
				  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=401; i<500; i++){
				  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
						  }else{
							  printf(" ");
						  }
			  }
			  printf("\n\n");
			  printf("\nquery_2\t\t\t\t\t");
			for(i=501;i<600; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=501;i<600; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=501;i<600; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=601;i<700; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=601;i<700; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=601;i<700; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=701;i<800; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=701;i<800; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=701;i<800; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=801;i<900; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=801;i<900; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=801;i<900; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='-' && aln->result_b[i]!='-' && aln->result_a[i]!='N'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=901; aln->result_b[i]!='\0'; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=901;aln->result_b[i]!='\0'; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=901;aln->result_b[i]!='\0'; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='-' && aln->result_b[i]!='-' && aln->result_a[i]!='N'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");*/
		//}
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		//printf("assigningscores to 2nd pair...\n");
		}
		FILE* site_scores_file;
		if ( print_all_nodes == 1 ){
			if (( site_scores_file = fopen("site_scores.txt","a")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
			//fprintf(site_scores_file,"%s\t",reverse_name);
		}
		// Initialize early termination state for this tree (second read of pair)
		type_of_PP best_score2 = -9999999999999999;
		int strikes2 = 0;
		type_of_PP strike_box_threshold2 = Cinterval * strike_box;
		type_of_PP pruning_threshold_calc2 = Cinterval * pruning_threshold;

		assignScores_Arr_paired(leaf_coordinates[match][0],rootArr[leaf_coordinates[match][0]],locQuery,positions,nodeScores,alength,match,print_all_nodes,site_scores_file,reverse_name,
		    early_termination, &best_score2, &strikes2, strike_box_threshold2, max_strikes,
		    enable_pruning, pruning_threshold_calc2);
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		//free(leaf_sequence);
		if ( print_all_nodes == 1 ){
			fclose(site_scores_file);
		}
		mm_allocator_free(wf_aligner->mm_allocator,pattern_alg);
		mm_allocator_free(wf_aligner->mm_allocator,ops_alg);
		mm_allocator_free(wf_aligner->mm_allocator,text_alg);
		}
		}
	}
	/* Delete the reused WFA aligner after all matches are done */
	wavefront_aligner_delete(wf_aligner);
	type_of_PP maximum=-9999999999999999;
	int minRoot=0;
	int minNode=0;
	int match_number=0;
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//printf("finding minimum score...\n");
	FILE* node_scores_file;
	if ( print_all_nodes == 1 ){
		if (( node_scores_file = fopen("scores_all_nodes.txt","w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
		fprintf(node_scores_file,"Tree_Number\tNode_Number\tScore\n");
	}
	for (i=0; i<number_of_matches;i++){
		for (j=leaf_coordinates[i][0]; j<leaf_coordinates[i][0]+1; j++){
			for(k=0; k<2*numspecArr[j]-1; k++){
				//printf("Tree %d Node %d Taxonomy (%s) Score: %lf\n",j,k,taxonomyArr[j][treeArr[j][k].taxIndex[0]][treeArr[j][k].taxIndex[1]],nodeScores[i][j][k]);
				if ( maximum < nodeScores[i][j][k]){
					maximum=nodeScores[i][j][k];
					match_number=i;
					minRoot=j;
					minNode=k;
				}
				if ( print_all_nodes == 1){
					fprintf(node_scores_file,"%d\t%d\t%lf\n",j,k,nodeScores[i][j][k]);
				}
			}
		}
	}
	if (print_all_nodes == 1){
		fclose(node_scores_file);
	}
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	//printf("match number: %d\n",match_number);
	//printf("minimum score: %lf\n",maximum);
	//printf("minimum root: %d\n",minRoot);
	//printf("minimum node: %d\n",minNode);
	//printf("C interval: %lf\n",Cinterval);
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//printf("clearing voteroot...\n");
	/*for(i=0; i<numberOfTotalRoots; i++){
		for(j=0;j<2*numspecArr[i]-1;j++){
			voteRoot[i][j]=0.0;
		}
	}*/
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	int index = 0;
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//printf("filling out voteroot...\n");
	for(i=0; i<number_of_matches; i++){
		//for(j=leaf_coordinates[i][0]; j<leaf_coordinates[i][0]+1; j++){
			for(k=0; k<2*numspecArr[leaf_coordinates[i][0]]-1; k++){
				if ( nodeScores[i][leaf_coordinates[i][0]][k] >= (maximum-Cinterval) && nodeScores[i][leaf_coordinates[i][0]][k] <= (maximum+Cinterval) ){
					//printf("Match : %d Min Root: %d Min node: %d, score: %lf\n",i,j,k,nodeScores[i][leaf_coordinates[i][0]][k]);
					voteRoot[leaf_coordinates[i][0]][k] = 1.0;
					index++;
				}
			}
		//}
	}
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));

	/* Best-leaf override: if the best-scoring leaf exceeds threshold and
	   total votes are low, override LCA by keeping only the best leaf's vote */
	if (best_leaf_max_votes > 0 && number_of_matches > 0) {
		type_of_PP bl_best = -9999999999999999;
		int bl_node = -1, bl_tree = -1;
		int bl_total = 0;
		for (i = 0; i < number_of_matches; i++) {
			int t = leaf_coordinates[i][0];
			if (t == -1) continue;
			int nn = 2*numspecArr[t]-1;
			for (k = 0; k < nn; k++) {
				if (treeArr[t][k].up[0] == -1 && treeArr[t][k].up[1] == -1) {
					if (nodeScores[i][t][k] > bl_best) {
						bl_best = nodeScores[i][t][k];
						bl_node = k; bl_tree = t;
					}
				}
				if (voteRoot[t][k] > 0.0) bl_total++;
			}
		}
		if (bl_node >= 0 && bl_best > best_leaf_threshold && bl_total < best_leaf_max_votes) {
			for (i = 0; i < number_of_matches; i++) {
				int t = leaf_coordinates[i][0];
				if (t == -1) continue;
				int nn = 2*numspecArr[t]-1;
				for (k = 0; k < nn; k++) voteRoot[t][k] = 0.0;
			}
			voteRoot[bl_tree][bl_node] = 1.0;
		}
	}

	/* Trace: report best-scoring leaf vs best-scoring node overall */
	if (g_trace_read[0] != '\0') {
		int is_traced = (g_trace_read[0] == '*' && g_trace_read[1] == '\0') ? 1 : strcmp(forward_name, g_trace_read) == 0;
		if (is_traced && number_of_matches > 0) {
			/* Find best leaf, second-best leaf, and count voted leaves/internal */
			type_of_PP best_leaf_score = -9999999999999999;
			type_of_PP second_leaf_score = -9999999999999999;
			int best_leaf_node = -1, best_leaf_tree = -1;
			int second_leaf_node = -1, second_leaf_tree = -1;
			int voted_leaves = 0, voted_internal = 0;
			for (i = 0; i < number_of_matches; i++) {
				int t = leaf_coordinates[i][0];
				int num_nodes = 2*numspecArr[t]-1;
				for (k = 0; k < num_nodes; k++) {
					int is_leaf = (treeArr[t][k].up[0] == -1 && treeArr[t][k].up[1] == -1);
					if (is_leaf) {
						type_of_PP s = nodeScores[i][t][k];
						if (s > best_leaf_score) {
							second_leaf_score = best_leaf_score;
							second_leaf_node = best_leaf_node;
							second_leaf_tree = best_leaf_tree;
							best_leaf_score = s;
							best_leaf_node = k;
							best_leaf_tree = t;
						} else if (s > second_leaf_score) {
							second_leaf_score = s;
							second_leaf_node = k;
							second_leaf_tree = t;
						}
					}
					if (voteRoot[t][k] > 0.0) {
						if (is_leaf) voted_leaves++;
						else voted_internal++;
					}
				}
			}
			int best_is_leaf = (minNode == best_leaf_node && minRoot == best_leaf_tree);
			fprintf(stderr, "  SCORING: max_node=tree%d:node%d score=%.2f is_leaf=%d\n",
				minRoot, minNode, (double)maximum, best_is_leaf);
			if (best_leaf_node >= 0) {
				int ti0 = treeArr[best_leaf_tree][best_leaf_node].taxIndex[0];
				int ti1 = treeArr[best_leaf_tree][best_leaf_node].taxIndex[1];
				fprintf(stderr, "  BEST_LEAF: tree%d:node%d score=%.2f gap=%.2f",
					best_leaf_tree, best_leaf_node, (double)best_leaf_score,
					(double)(maximum - best_leaf_score));
				if (ti0 != -1 && ti1 != -1)
					fprintf(stderr, " tax=%s", taxonomyArr[best_leaf_tree][ti0][ti1]);
				fprintf(stderr, "\n");
			}
			if (second_leaf_node >= 0) {
				int ti0 = treeArr[second_leaf_tree][second_leaf_node].taxIndex[0];
				int ti1 = treeArr[second_leaf_tree][second_leaf_node].taxIndex[1];
				fprintf(stderr, "  2ND_LEAF: tree%d:node%d score=%.2f gap=%.2f",
					second_leaf_tree, second_leaf_node, (double)second_leaf_score,
					(double)(best_leaf_score - second_leaf_score));
				if (ti0 != -1 && ti1 != -1)
					fprintf(stderr, " tax=%s", taxonomyArr[second_leaf_tree][ti0][ti1]);
				fprintf(stderr, "\n");
			}
			fprintf(stderr, "  VOTES: %d leaves + %d internal = %d total (Cinterval=%.1f)\n",
				voted_leaves, voted_internal, voted_leaves + voted_internal, (double)Cinterval);
		}
	}

	minimum_score[0] = maximum;
	minimum_score[1] = forward_mismatch;
	minimum_score[2] = reverse_mismatch;
	//printf("%lf\t",minimum);
	//free(leaf_sequence);
	//free(positionsInRoot);
}
void place_paired_with_nw( char *query_1, char *query_2, char **rootSeqs, int numberOfTotalRoots, int *positions, char *locQuery, nw_aligner_t *nw, alignment_t *aln, scoring_t *scoring, type_of_PP ***nodeScores, double **voteRoot, int number_of_matches , int **leaf_coordinates, int paired, type_of_PP* minimum_score, char *alignments_dir, char *forward_name, char *reverse_name, int print_alignments, char *leaf_sequence, int *positionsInRoot, int maxNumSpec, int* starts_forward, char** cigars_forward, int* starts_reverse, char** cigars_reverse, int print_alignments_to_file, int use_leaf_portion, int padding, int max_query_length, int max_numbase, int print_all_nodes, int early_termination, type_of_PP strike_box, int max_strikes, int enable_pruning, type_of_PP pruning_threshold, type_of_PP best_leaf_threshold, int best_leaf_max_votes){
	int i, j, k, node, match;
	type_of_PP forward_mismatch, reverse_mismatch;
	forward_mismatch=0;
	reverse_mismatch=0;
	for(match=0; match<number_of_matches; match++){
		{
			int clear_len = use_leaf_portion ? (max_query_length+max_query_length+2*padding+1) : (max_query_length+max_numbase+1);
			memset(positions, -1, clear_len * sizeof(int));
			memset(locQuery, '\0', clear_len);
			memset(leaf_sequence, '\0', clear_len);
			memset(positionsInRoot, -1, clear_len * sizeof(int));
		}
		int query_length = strlen(query_1);
		if (leaf_coordinates[match][0] == -1 || leaf_coordinates[match][1] == -1) {
			continue;  /* Skip invalid leaf coordinates */
		}
		if (use_leaf_portion==1){
			if ( cigars_forward[match][0] == '*'){ break; }
			if ( cigars_forward[match][0] == '\0'){ break; }
		}
		if (use_leaf_portion == 1 && starts_forward[match] != -1){
			int start_position = getStartPosition(starts_forward[match],leaf_coordinates[match][0],leaf_coordinates[match][1],padding);
			int end_position = getEndPosition(cigars_forward[match],leaf_coordinates[match][0],leaf_coordinates[match][1],start_position+padding,padding);
			getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,start_position,end_position);
		}else{
			getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,0,numbaseArr[leaf_coordinates[match][0]]);
		}
		int leaf_length = strlen(leaf_sequence);
		if (leaf_length > 0){
		/*printf("leaf sequence: ");
		for(i=0; i<strlen(leaf_sequence); i++){
			printf("%c",leaf_sequence[i]);
		}*/
		//printf("\n");
		//printf("aligning...\n");
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		needleman_wunsch_align(leaf_sequence, query_1, scoring, nw, aln);
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		int alength=0;
		int keepTrackOfPosInRoot=0;
		int loop=0;
		int lengthOfResultA=0;
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (aln->result_a[i]=='-'){*/
				/*positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;*/
				//keepTrackOfPosInRoot++;
			/*}else if (aln->result_a[i]=='N'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;
				keepTrackOfPosInRoot++;
			}else{
				if (alength != 0 && aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
				//if ( aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
					locQuery[alength]=aln->result_b[loop];
					positions[alength]=keepTrackOfPosInRoot;
					alength++;
				}
				if ( aln->result_b[loop] == '-' && aln->result_a[loop] != '-'){
					keepTrackOfPosInRoot++;
					loop++;
				}else{
					if (aln->result_a[i]!='-'){
						locQuery[alength]=aln->result_b[loop];
						positions[alength]=keepTrackOfPosInRoot;
						keepTrackOfPosInRoot++;
						alength++;
						loop++;
						lengthOfResultA++;
					}
				}
			}
		}*/
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (keepTrackOfPosInRoot > leaf_length){
				positions[0]=-1;
				break;
			}
			//if (aln->result_a[i] != '-' && aln->result_b[i] != '-'){
			if (aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				if (aln->result_a[i] != aln->result_b[i] && match==0){
					forward_mismatch++;
				}
			}
			if (alength != 0 && lengthOfResultA<query_length && aln->result_b[i] =='-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				keepTrackOfPosInRoot++;
				if (aln->result_a[i] != aln->result_b[i] && match==0 && aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
					forward_mismatch++;
				}
			}else if (aln->result_a[i] != '-'){
			//}else if (aln->result_a[i] != 'N'){
				keepTrackOfPosInRoot++;
			}
		}
		locQuery[alength]='\0';*/
		int counterInRoot=0;
		int counterinB=0;
			for(i=0; aln->result_a[i] != '\0'; i++){
				if (aln->result_a[i] != '-' && aln->result_b[i] != '-'){
					positions[alength] = positionsInRoot[counterInRoot];
					locQuery[alength] = aln->result_b[i];
					alength++;
					if (aln->result_a[i] != aln->result_b[i] && match==0){
						forward_mismatch++;
					}
				}
				if (aln->result_a[i] == '-' && aln->result_b[i] != '-' && alength>0 && counterinB<query_length){
					positions[alength]=-1;
					locQuery[alength] = aln->result_b[i];
					alength++;	
				}
				if (aln->result_a[i] != '-' && aln->result_b[i] == '-' && alength > 0 && counterinB<query_length){
					positions[alength]= positionsInRoot[counterInRoot];
					locQuery[alength] = aln->result_b[i];
					alength++;
				}
				if (aln->result_a[i] != '-'){
					counterInRoot++;
				}
				if (aln->result_b[i] != '-'){
					counterinB++;
				}
			}
		locQuery[alength]='\0';
		//if (access(alignmentFileName, F_OK ) != -1 && match==0){
		//	printToFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, forward_name);
		//}else if (match==0){
		//	createNewFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, forward_name);
		//}
		if (match==0 && print_alignments_to_file==1){
			char alignmentFileName[1000];
			snprintf(alignmentFileName,1000,"%s/%s.fasta",alignments_dir,treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			//snprintf(alignmentFileName,1000,"%s",alignments_dir);
			if ( access(alignmentFileName, F_OK ) != -1 ){
				printToFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, forward_name, query_length, positionsInRoot);
			}else{
				createNewFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, forward_name, query_length, positionsInRoot);
			}
		}
		if (print_alignments==1){
			printf("Using Needleman-Wunsch\n");
			printf("\n");
			printf("%s\n",forward_name);
			  printf("query_1\t\t\t\t\t");
			  for(i=0; i<100; i++){
			  		printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=0; i<100; i++){
					printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=0; i<100; i++){
			  		if ( aln->result_a[i] == aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  			printf("|");
			  		}else if (aln->result_a[i] != aln->result_b[i] && aln->result_a[i] !='N' && aln->result_b[i] != '-'){
			  			printf("*");
			  		}else{
			  			printf(" ");
			  		}
			  }
			  printf("\n\n");
			printf("%s\n",forward_name);
			  printf("query_1\t\t\t\t\t");
			  for(i=100;i<200;i++){
				  	printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=100;i<200;i++){
					printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=100;i<200;i++){
					  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
						  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
					  }else{
						  printf(" ");
					  }
			  }
			  printf("\n\n");
			printf("%s\n",forward_name);
			  printf("query_1\t\t\t\t\t");
			  for(i=200; i<300; i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=200;i<300;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=200;i<300;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
			  printf("*");
			  }else{
			  printf(" ");
			  }
			  }
			  printf("\n\n");
			printf("%s\n",forward_name);
			  printf("query_1\t\t\t\t\t");
			  for(i=300; i<400 ; i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=300;i<400; i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=300;i<400;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
				printf("*");
			  }else{
				  printf(" ");
			  }
			  }
			  printf("\n\n");
			  /*
			  printf("query_1\t\t\t\t\t");
			  for(i=401; i<500 ; i++){
				printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=401; i<500 ; i++){
				  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=401; i<500 ; i++){
				  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
						  }else{
							  printf(" ");
						  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=501; i<600; i++){
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=501;i<600; i++){
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=501;i<600;i++){
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
				  printf("*");
			  }else{
				  printf(" ");
			  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=601; i<700; i++){
				  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=601;i<700; i++){
				  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=601;i<700;i++){
				  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
						  }else{
							  printf(" ");
						  }
			  }
			  printf("\n\n");
			printf("query_1\t\t\t\t\t");
			for(i=701; i<800; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=701;i<800; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=701;i<800;i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
					}else{
						printf(" ");
					}
			}
			printf("\n\n");
			printf("query_1\t\t\t\t\t");
			for(i=801;i<900; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=801;i<900; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=801;i<900;i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					 printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						   printf("*");
					  }else{
						  printf(" ");
					  }
			}
			 printf("\n\n");
			printf("query_1\t\t\t\t\t");
			for(i=901; aln->result_b[i] != '\0'; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=901;aln->result_b[i]!='\0'; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=901;aln->result_b[i]!='\0';i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
					}else{
						printf(" ");
					}
			}
			printf("\n\n");*/
		}
		//assignScores_Arr_paired(topRoots[rootNum],rootArr[topRoots[rootNum]],locQuery, positions, nodeScores, alength);
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		//printf("assignScores_Arr_paired...\n");
		FILE* site_scores_file;
		if ( print_all_nodes == 1 && match > 0 && access("site_scores.txt", F_OK ) != -1 ){
			if (( site_scores_file = fopen("site_scores.txt","a")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
			//fprintf(site_scores_file,"%s\t",forward_name);
		} else if ( print_all_nodes == 1 && match ==0 ){
			if (( site_scores_file = fopen("site_scores.txt","w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
			fprintf(site_scores_file,"Readname\tTree_Number\tNode_Number\tSite\tScore\n");
			//fprintf(site_scores_file,"%s\t",forward_name);
		}
		// Initialize early termination state for this tree
		type_of_PP best_score = -9999999999999999;
		int strikes = 0;
		type_of_PP strike_box_threshold = Cinterval * strike_box;
		type_of_PP pruning_threshold_calc = Cinterval * pruning_threshold;

		assignScores_Arr_paired(leaf_coordinates[match][0],rootArr[leaf_coordinates[match][0]],locQuery, positions, nodeScores, alength, match, print_all_nodes, site_scores_file,forward_name,
		    early_termination, &best_score, &strikes, strike_box_threshold, max_strikes,
		    enable_pruning, pruning_threshold_calc);
		if ( print_all_nodes == 1){
			fclose(site_scores_file);
		}
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	}
	}
	if (paired==1){
	for(match=0; match<number_of_matches; match++){
		{
			int clear_len = use_leaf_portion ? (max_query_length+2*padding+1) : (max_query_length+max_numbase+1);
			memset(positions, -1, clear_len * sizeof(int));
			memset(locQuery, '\0', clear_len);
			memset(leaf_sequence, '\0', clear_len);
			memset(positionsInRoot, -1, clear_len * sizeof(int));
		}
		int query_length = strlen(query_2);
		if (leaf_coordinates[match][0] == -1 || leaf_coordinates[match][1] == -1) {
			continue;  /* Skip invalid leaf coordinates */
		}
		//needleman_wunsch_align(rootSeqs[topRoots[rootNum]], query_2, scoring, nw, aln);
		//char *leaf_sequence = (char *)malloc(numbaseArr[leaf_coordinates[rootNum][0]]*sizeof(char));
		//getSequenceInNode(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence);
		if (use_leaf_portion==1){
			if ( cigars_reverse[match][0] == '*'){
				break;
			}
			if ( cigars_forward[match][0] == '\0'){ break; }
		}
		if (use_leaf_portion == 1 ){
			if ( starts_reverse[match]!=1){
				int start_position = getStartPosition(starts_reverse[match],leaf_coordinates[match][0],leaf_coordinates[match][1],padding);
				int end_position = getEndPosition(cigars_reverse[match],leaf_coordinates[match][0],leaf_coordinates[match][1],start_position+padding,padding);
				getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,start_position,end_position);
			}
		}else{
			getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,0,numbaseArr[leaf_coordinates[match][0]]);
		}
		int leaf_length = strlen(leaf_sequence);
		if(leaf_length > 0){
		//printf("leaf sequence: ");
		//for(i=0; i<strlen(leaf_sequence); i++){
		//	printf("%c",leaf_sequence[i]);
		//}
		//printf("\n");
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		//printf("aligning 2nd pair...\n");
		needleman_wunsch_align(leaf_sequence, query_2, scoring, nw, aln);
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		/*printf("leaf sequence: ");
		for(i=0; i<strlen(leaf_sequence); i++){
			printf("%c",leaf_sequence[i]);
		}*/
		int alength=0;
		int keepTrackOfPosInRoot=0;
		int loop=0;
		int lengthOfResultA=0;
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (aln->result_a[i]=='-'){*/
				/*positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;*/
				//keepTrackOfPosInRoot++;
			/*}else if (aln->result_a[i]=='N'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;
				keepTrackOfPosInRoot++;
			}else{
				if (alength != 0 && aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
				//if ( aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
					locQuery[alength]=aln->result_b[loop];
					positions[alength]=keepTrackOfPosInRoot;
					alength++;
				}
				if ( aln->result_b[loop] == '-' && aln->result_a[loop] != '-'){
					keepTrackOfPosInRoot++;
					loop++;
				}else{
					if (aln->result_a[i]!='-'){
						locQuery[alength]=aln->result_b[loop];
						positions[alength]=keepTrackOfPosInRoot;
						keepTrackOfPosInRoot++;
						alength++;
						loop++;
						lengthOfResultA++;
					}
				}
			}
		}*/
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (keepTrackOfPosInRoot > leaf_length){
				positions[0]=-1;
				break;
			}
			//if (aln->result_a[i] != '-' && aln->result_b[i] != '-'){
			if (aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				if (aln->result_a[i] != aln->result_b[i] && match==0){
					reverse_mismatch++;
				}
			}
			if (alength != 0 && lengthOfResultA<query_length && aln->result_b[i]=='-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				keepTrackOfPosInRoot++;
				if (aln->result_a[i] != aln->result_b[i] && match==0 && aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
					reverse_mismatch++;
				}
			}else if (aln->result_a[i] != '-'){
			//}else if (aln->result_a[i] != 'N'){
				keepTrackOfPosInRoot++;
			}
		}
		locQuery[alength]='\0';*/
		int counterInRoot=0;
		int counterinB=0;
			for(i=0; aln->result_a[i] != '\0'; i++){
				if (aln->result_a[i] != '-' && aln->result_b[i] != '-'){
					positions[alength] = positionsInRoot[counterInRoot];
					locQuery[alength] = aln->result_b[i];
					alength++;
					if (aln->result_a[i] != aln->result_b[i] && match==0){
						reverse_mismatch++;
					}
				}
				if (aln->result_a[i] == '-' && aln->result_b[i] != '-' && alength>0 && counterinB<query_length){
					positions[alength]=-1;
					locQuery[alength] = aln->result_b[i];
					alength++;
				}
				if (aln->result_a[i] != '-' && aln->result_b[i] == '-' && alength > 0 && counterinB <query_length){
					positions[alength]= positionsInRoot[counterInRoot];
					locQuery[alength] = aln->result_b[i];
					alength++;
				}
				if (aln->result_a[i] != '-'){
					counterInRoot++;
				}
				if (aln->result_b[i] != '-'){
					counterinB++;
				}
			}
		locQuery[alength]='\0';
		//if (access(alignmentFileName, F_OK ) != -1 && match==0){
		//	printToFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, reverse);
		//}else if (match==0){
		//	createNewFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, readname);
		//}
		if (match==0 && print_alignments_to_file==1){
			char alignmentFileName[1000];
			snprintf(alignmentFileName,1000,"%s/%s.fasta",alignments_dir,treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			if ( access(alignmentFileName, F_OK ) != -1 ){
				printToFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, reverse_name, query_length, positionsInRoot);
			}else{
				createNewFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, reverse_name, query_length, positionsInRoot);
			}
		}
		if (print_alignments==1){
			printf("\n");
			printf("%s\n",reverse_name);
			printf("query_2\t\t\t\t\t");
			for(i=0; i<100; i++){
					printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=0; i<100; i++){
					printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t");
			for(i=0; i<100; i++){
					if ( aln->result_a[i] == aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
						printf("|");
					}else if (aln->result_a[i] !='N' && aln->result_b[i] != '-'){
						printf("*");
					}else{
						printf(" ");
					}
			}
			printf("\n\n");
			printf("%s\n",reverse_name);
			printf("query_2\t\t\t\t\t");
			  for(i=101;i<200;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=101;i<200;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=101;i<200;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
			  printf("*");
			  }else{
			  printf(" ");
			  }
			  }
			  printf("\n\n");
			printf("%s\n",reverse_name);
			  printf("query_2\t\t\t\t\t");
			  for(i=201; i<300; i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=201;i<300;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=201;i<300;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
			  printf("*");
			  }else{
			  printf(" ");
			  }
			  }
			  printf("\n\n");
			printf("%s\n",reverse_name);
			  printf("\nquery_2\t\t\t\t\t");
			  for(i=301; i<400; i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=301;i<400; i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=301;i<400;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
				  printf("*");
			  }else{
				  printf(" ");
			  }
			  }
			  printf("\n\n");
			  /*
			  printf("\nquery_2\t\t\t\t\t");
			  for(i=401; i<500; i++){
				  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=401; i<500; i++){
				  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=401; i<500; i++){
				  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
						  }else{
							  printf(" ");
						  }
			  }
			  printf("\n\n");
			  printf("\nquery_2\t\t\t\t\t");
			for(i=501;i<600; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=501;i<600; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=501;i<600; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=601;i<700; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=601;i<700; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=601;i<700; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=701;i<800; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=701;i<800; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=701;i<800; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=801;i<900; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=801;i<900; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=801;i<900; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='-' && aln->result_b[i]!='-' && aln->result_a[i]!='N'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=901; aln->result_b[i]!='\0'; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=901;aln->result_b[i]!='\0'; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=901;aln->result_b[i]!='\0'; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='-' && aln->result_b[i]!='-' && aln->result_a[i]!='N'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");*/
		//}
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		//printf("assigningscores to 2nd pair...\n");
		}
		FILE* site_scores_file;
		if ( print_all_nodes == 1 ){
			if (( site_scores_file = fopen("site_scores.txt","a")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
			//fprintf(site_scores_file,"%s\t",reverse_name);
		}
		// Initialize early termination state for this tree (second read of pair)
		type_of_PP best_score2 = -9999999999999999;
		int strikes2 = 0;
		type_of_PP strike_box_threshold2 = Cinterval * strike_box;
		type_of_PP pruning_threshold_calc2 = Cinterval * pruning_threshold;

		assignScores_Arr_paired(leaf_coordinates[match][0],rootArr[leaf_coordinates[match][0]],locQuery,positions,nodeScores,alength,match,print_all_nodes,site_scores_file,reverse_name,
		    early_termination, &best_score2, &strikes2, strike_box_threshold2, max_strikes,
		    enable_pruning, pruning_threshold_calc2);
		if ( print_all_nodes == 1){
			fclose(site_scores_file);
		}
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		//free(leaf_sequence);
	}
		}
	}
	type_of_PP maximum=-9999999999999999;
	int minRoot=0;
	int minNode=0;
	int match_number=0;
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//printf("finding minimum score...\n");
	FILE* node_scores_file;
	if ( print_all_nodes == 1 ){
		if (( node_scores_file = fopen("scores_all_nodes.txt","w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
		fprintf(node_scores_file,"Tree_Number\tNode_Number\tScore\n");
	}
	for (i=0; i<number_of_matches;i++){
		for (j=leaf_coordinates[i][0]; j<leaf_coordinates[i][0]+1; j++){
			for(k=0; k<2*numspecArr[j]-1; k++){
				if ( maximum < nodeScores[i][j][k]){
					maximum=nodeScores[i][j][k];
					match_number=i;
					minRoot=j;
					minNode=k;
				}
				if ( print_all_nodes == 1){
					fprintf(node_scores_file,"%d\t%d\t%lf\n",j,k,nodeScores[i][j][k]);
				}
			}
		}
	}
	if (print_all_nodes == 1){
		fclose(node_scores_file);
	}
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	//printf("match number: %d\n",match_number);
	//printf("minimum score: %lf\n",maximum);
	//printf("minimum root: %d\n",minRoot);
	//printf("minimum node: %d\n",minNode);
	//printf("C interval: %lf\n",Cinterval);
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//printf("clearing voteroot...\n");
	/*for(i=0; i<numberOfTotalRoots; i++){
		for(j=0;j<2*numspecArr[i]-1;j++){
			voteRoot[i][j]=0.0;
		}
	}*/
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	int index = 0;
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//printf("filling out voteroot...\n");
	for(i=0; i<number_of_matches; i++){
		//for(j=leaf_coordinates[i][0]; j<leaf_coordinates[i][0]+1; j++){
			for(k=0; k<2*numspecArr[leaf_coordinates[i][0]]-1; k++){
				if ( nodeScores[i][leaf_coordinates[i][0]][k] >= (maximum-Cinterval) && nodeScores[i][leaf_coordinates[i][0]][k] <= (maximum+Cinterval) ){
					//printf("Match : %d Min Root: %d Min node: %d, score: %lf\n",i,j,k,nodeScores[i][leaf_coordinates[i][0]][k]);
					voteRoot[leaf_coordinates[i][0]][k] = 1.0;
					index++;
				}
			}
		//}
	}
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));

	/* Best-leaf override (NW path) */
	if (best_leaf_max_votes > 0 && number_of_matches > 0) {
		type_of_PP bl_best = -9999999999999999;
		int bl_node = -1, bl_tree = -1;
		int bl_total = 0;
		for (i = 0; i < number_of_matches; i++) {
			int t = leaf_coordinates[i][0];
			if (t == -1) continue;
			int nn = 2*numspecArr[t]-1;
			for (k = 0; k < nn; k++) {
				if (treeArr[t][k].up[0] == -1 && treeArr[t][k].up[1] == -1) {
					if (nodeScores[i][t][k] > bl_best) {
						bl_best = nodeScores[i][t][k];
						bl_node = k; bl_tree = t;
					}
				}
				if (voteRoot[t][k] > 0.0) bl_total++;
			}
		}
		if (bl_node >= 0 && bl_best > best_leaf_threshold && bl_total < best_leaf_max_votes) {
			for (i = 0; i < number_of_matches; i++) {
				int t = leaf_coordinates[i][0];
				if (t == -1) continue;
				int nn = 2*numspecArr[t]-1;
				for (k = 0; k < nn; k++) voteRoot[t][k] = 0.0;
			}
			voteRoot[bl_tree][bl_node] = 1.0;
		}
	}

	/* Trace: report best-scoring leaf vs best-scoring node overall (NW path) */
	if (g_trace_read[0] != '\0') {
		int is_traced_nw = (g_trace_read[0] == '*' && g_trace_read[1] == '\0') ? 1 : strcmp(forward_name, g_trace_read) == 0;
		if (is_traced_nw && number_of_matches > 0) {
			type_of_PP best_leaf_score_nw = -9999999999999999;
			type_of_PP second_leaf_score_nw = -9999999999999999;
			int best_leaf_node_nw = -1, best_leaf_tree_nw = -1;
			int second_leaf_node_nw = -1, second_leaf_tree_nw = -1;
			int voted_leaves_nw = 0, voted_internal_nw = 0;
			for (i = 0; i < number_of_matches; i++) {
				int t = leaf_coordinates[i][0];
				if (t == -1) continue;
				int num_nodes = 2*numspecArr[t]-1;
				for (k = 0; k < num_nodes; k++) {
					int is_lf = (treeArr[t][k].up[0] == -1 && treeArr[t][k].up[1] == -1);
					if (is_lf) {
						type_of_PP s = nodeScores[i][t][k];
						if (s > best_leaf_score_nw) {
							second_leaf_score_nw = best_leaf_score_nw;
							second_leaf_node_nw = best_leaf_node_nw;
							second_leaf_tree_nw = best_leaf_tree_nw;
							best_leaf_score_nw = s;
							best_leaf_node_nw = k;
							best_leaf_tree_nw = t;
						} else if (s > second_leaf_score_nw) {
							second_leaf_score_nw = s;
							second_leaf_node_nw = k;
							second_leaf_tree_nw = t;
						}
					}
					if (voteRoot[t][k] > 0.0) {
						if (is_lf) voted_leaves_nw++;
						else voted_internal_nw++;
					}
				}
			}
			int best_is_leaf_nw = (minNode == best_leaf_node_nw && minRoot == best_leaf_tree_nw);
			fprintf(stderr, "  SCORING: max_node=tree%d:node%d score=%.2f is_leaf=%d\n",
				minRoot, minNode, (double)maximum, best_is_leaf_nw);
			if (best_leaf_node_nw >= 0) {
				int ti0 = treeArr[best_leaf_tree_nw][best_leaf_node_nw].taxIndex[0];
				int ti1 = treeArr[best_leaf_tree_nw][best_leaf_node_nw].taxIndex[1];
				fprintf(stderr, "  BEST_LEAF: tree%d:node%d score=%.2f gap=%.2f",
					best_leaf_tree_nw, best_leaf_node_nw, (double)best_leaf_score_nw,
					(double)(maximum - best_leaf_score_nw));
				if (ti0 != -1 && ti1 != -1)
					fprintf(stderr, " tax=%s", taxonomyArr[best_leaf_tree_nw][ti0][ti1]);
				fprintf(stderr, "\n");
			}
			if (second_leaf_node_nw >= 0) {
				int ti0 = treeArr[second_leaf_tree_nw][second_leaf_node_nw].taxIndex[0];
				int ti1 = treeArr[second_leaf_tree_nw][second_leaf_node_nw].taxIndex[1];
				fprintf(stderr, "  2ND_LEAF: tree%d:node%d score=%.2f gap=%.2f",
					second_leaf_tree_nw, second_leaf_node_nw, (double)second_leaf_score_nw,
					(double)(best_leaf_score_nw - second_leaf_score_nw));
				if (ti0 != -1 && ti1 != -1)
					fprintf(stderr, " tax=%s", taxonomyArr[second_leaf_tree_nw][ti0][ti1]);
				fprintf(stderr, "\n");
			}
			fprintf(stderr, "  VOTES: %d leaves + %d internal = %d total (Cinterval=%.1f)\n",
				voted_leaves_nw, voted_internal_nw, voted_leaves_nw + voted_internal_nw, (double)Cinterval);
		}
	}

	minimum_score[0] = maximum;
	minimum_score[1] = forward_mismatch;
	minimum_score[2] = reverse_mismatch;
	//printf("%lf\t",minimum);
	//free(leaf_sequence);
	//free(positionsInRoot);
}
