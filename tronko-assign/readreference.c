#include "readreference.h"
#include "crash_debug.h"
#include "logger.h"
#include "resource_monitor.h"
#include <stdint.h>  // For uint8_t, uint16_t, uint32_t, uint64_t in format detection

// External TSV log file handle (set in main)
extern FILE *g_tsv_log_file;

// Little-endian read helpers for binary format
// (Matches tronko-convert/utils.c naming convention)
static uint16_t read_u16(FILE *fp) {
    uint8_t b[2];
    if (fread(b, 1, 2, fp) != 2) {
        LOG_ERROR("Unexpected end of file reading uint16");
        return 0;
    }
    return (uint16_t)b[0] | ((uint16_t)b[1] << 8);
}

static uint32_t read_u32(FILE *fp) {
    uint8_t b[4];
    if (fread(b, 1, 4, fp) != 4) {
        LOG_ERROR("Unexpected end of file reading uint32");
        return 0;
    }
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

static int32_t read_i32(FILE *fp) {
    return (int32_t)read_u32(fp);
}

static uint64_t read_u64(FILE *fp) {
    uint8_t b[8];
    if (fread(b, 1, 8, fp) != 8) {
        LOG_ERROR("Unexpected end of file reading uint64");
        return 0;
    }
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) {
        v |= ((uint64_t)b[i] << (i * 8));
    }
    return v;
}

// Gzip-compatible little-endian read helpers for binary format
static uint16_t gz_read_u16(gzFile gz) {
    uint8_t b[2];
    if (gzread(gz, b, 2) != 2) {
        LOG_ERROR("Unexpected end of file reading uint16 (gzip)");
        return 0;
    }
    return (uint16_t)b[0] | ((uint16_t)b[1] << 8);
}

static uint32_t gz_read_u32(gzFile gz) {
    uint8_t b[4];
    if (gzread(gz, b, 4) != 4) {
        LOG_ERROR("Unexpected end of file reading uint32 (gzip)");
        return 0;
    }
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

static int32_t gz_read_i32(gzFile gz) {
    return (int32_t)gz_read_u32(gz);
}

static uint64_t gz_read_u64(gzFile gz) {
    uint8_t b[8];
    if (gzread(gz, b, 8) != 8) {
        LOG_ERROR("Unexpected end of file reading uint64 (gzip)");
        return 0;
    }
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) {
        v |= ((uint64_t)b[i] << (i * 8));
    }
    return v;
}

/**
 * Detect reference file format (binary or text, optionally gzipped)
 * Returns: FORMAT_BINARY, FORMAT_BINARY_GZIPPED, FORMAT_TEXT, or FORMAT_UNKNOWN
 */
int detect_reference_format(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        LOG_ERROR("Cannot open file for format detection: %s", filename);
        return FORMAT_UNKNOWN;
    }

    uint8_t magic[4];
    size_t n = fread(magic, 1, 4, fp);
    fclose(fp);

    if (n < 4) {
        LOG_ERROR("File too small for format detection: %s", filename);
        return FORMAT_UNKNOWN;
    }

    // Check for binary format magic: 0x89 'T' 'R' 'K'
    if (magic[0] == TRONKO_MAGIC_0 && magic[1] == TRONKO_MAGIC_1 &&
        magic[2] == TRONKO_MAGIC_2 && magic[3] == TRONKO_MAGIC_3) {
        LOG_DEBUG("Detected binary format: %s", filename);
        return FORMAT_BINARY;
    }

    // Check for gzip magic: 0x1f 0x8b
    if (magic[0] == 0x1f && magic[1] == 0x8b) {
        // Gzipped file - decompress first 4 bytes to check content type
        gzFile gz = gzopen(filename, "rb");
        if (gz) {
            uint8_t inner_magic[4];
            int bytes_read = gzread(gz, inner_magic, 4);
            gzclose(gz);

            if (bytes_read == 4) {
                if (inner_magic[0] == TRONKO_MAGIC_0 && inner_magic[1] == TRONKO_MAGIC_1 &&
                    inner_magic[2] == TRONKO_MAGIC_2 && inner_magic[3] == TRONKO_MAGIC_3) {
                    LOG_DEBUG("Detected gzipped binary format: %s", filename);
                    return FORMAT_BINARY_GZIPPED;
                }
            }
        }
        LOG_DEBUG("Detected gzipped text format: %s", filename);
        return FORMAT_TEXT;  // Gzipped text
    }

    // Check for zstd magic: 0x28 0xb5 0x2f 0xfd
    if (magic[0] == 0x28 && magic[1] == 0xb5 && magic[2] == 0x2f && magic[3] == 0xfd) {
        // Zstd compressed - use streaming decompression to check for TRKB magic
        FILE *zfp = fopen(filename, "rb");
        if (zfp) {
            ZSTD_DCtx *dctx = ZSTD_createDCtx();
            if (dctx) {
                size_t in_buf_size = ZSTD_DStreamInSize();
                size_t out_buf_size = ZSTD_DStreamOutSize();
                uint8_t *in_buf = malloc(in_buf_size);
                uint8_t *out_buf = malloc(out_buf_size);

                if (in_buf && out_buf) {
                    size_t n_read = fread(in_buf, 1, in_buf_size, zfp);
                    if (n_read > 0) {
                        ZSTD_inBuffer in = { in_buf, n_read, 0 };
                        ZSTD_outBuffer out = { out_buf, out_buf_size, 0 };

                        size_t ret = ZSTD_decompressStream(dctx, &out, &in);
                        if (!ZSTD_isError(ret) && out.pos >= 4) {
                            if (out_buf[0] == TRONKO_MAGIC_0 && out_buf[1] == TRONKO_MAGIC_1 &&
                                out_buf[2] == TRONKO_MAGIC_2 && out_buf[3] == TRONKO_MAGIC_3) {
                                LOG_DEBUG("Detected zstd-compressed binary format: %s", filename);
                                free(in_buf);
                                free(out_buf);
                                ZSTD_freeDCtx(dctx);
                                fclose(zfp);
                                return FORMAT_BINARY_ZSTD;
                            }
                        }
                    }
                }

                if (in_buf) free(in_buf);
                if (out_buf) free(out_buf);
                ZSTD_freeDCtx(dctx);
            }
            fclose(zfp);
        }
        LOG_WARN("Zstd file does not contain TRKB data: %s", filename);
        return FORMAT_UNKNOWN;
    }

    // Assume text if first byte is ASCII digit (numberOfTrees line)
    if (magic[0] >= '0' && magic[0] <= '9') {
        LOG_DEBUG("Detected plain text format: %s", filename);
        return FORMAT_TEXT;
    }

    LOG_WARN("Unknown file format: %s", filename);
    return FORMAT_UNKNOWN;
}
int readInXNumberOfLines_fastq(int numberOfLinesToRead, gzFile query_reads, int whichPair, Options opt, int max_query_length, int max_readname_length,int first_iter){
	char* buffer;
	char* query;
	char* reverse;
	int buffer_size = 0;
	if ( max_query_length > max_readname_length ){
		buffer_size = max_query_length + 2;
		buffer = (char *)malloc(sizeof(char)*(max_query_length+2));
	}else{
		buffer_size = max_readname_length + 2;
		buffer = (char *)malloc(sizeof(char)*(max_readname_length+2));
	}
	char* s;
	char seqname[max_readname_length+1];
	int size=0;
	int i=0;
	int iter=0;
	int next=0;
	query = (char *)malloc(sizeof(char)*max_query_length+2);
	reverse = (char *)malloc(sizeof(char)*max_query_length+2);
	for(i=0; i<max_query_length+2; i++){
		query[i] = '\0';
		reverse[i] = '\0';
	}
	int line_number_fastq = 0;
	// Set the current file being processed based on whichPair
	const char* current_filename = NULL;
	if (whichPair == 1) {
		current_filename = opt.read1_file;
	} else if (whichPair == 2) {
		current_filename = opt.read2_file;
	}
	if (current_filename) {
		crash_set_current_file(current_filename);
	}
	
	while(gzgets(query_reads,buffer,buffer_size)!=NULL){
		line_number_fastq++;
		crash_set_current_file_line(current_filename, line_number_fastq);  // Update current file and line number
		s = strtok(buffer,"\n");
		size = strlen(s);
		if ( first_iter == 1 ){
			if ( buffer[0] != '@' ){
				printf("Query reads are not in FASTQ format\n");
				exit(-1);
			}
			first_iter=0;
		}
		
		// Check for data corruption patterns in FASTQ
		if (buffer[0] == '@') {
			// Header line - validate format
			if (size < 2) {
				crash_flag_corruption(current_filename, line_number_fastq, "Empty FASTQ header");
			}
		} else if (buffer[0] != '+' && buffer[0] != '>') {
			// Sequence or quality line - check for corruption patterns
			if (buffer[0] == ' ' || buffer[0] == '\t') {
				crash_flag_corruption(current_filename, line_number_fastq, "FASTQ line starts with whitespace");
			}
			if (size > 0 && size < 10) {  // Very short sequences are suspicious
				crash_flag_corruption(current_filename, line_number_fastq, "Suspiciously short FASTQ sequence");
			}
		}
		if ( buffer[0] == '@' && whichPair==1){
			for(i=1; i<size; i++){
				if (buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			memset(pairedQueryMat->forward_name[iter],'\0',max_readname_length);
			strcpy(pairedQueryMat->forward_name[iter],seqname);
			memset(seqname,'\0',max_readname_length);
			next=1;
		}else if ( buffer[0] == '@' && whichPair==2){
			for(i=1; i<size; i++){
				if(buffer[i]==' '){ buffer[i]='_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			char tempname[max_readname_length];
			memset(tempname,'\0',max_readname_length);
			for(i=0; i<size-1; i++){
				if ( pairedQueryMat->forward_name[iter][i] == '1' && pairedQueryMat->forward_name[iter][i-1] == '_'){
					tempname[i] = '2';
				}else{
					tempname[i] = seqname[i];
				}
			}
			tempname[size-1]='\0';
			int skipped=iter;
			/*while(strcmp(tempname,seqname)!=0){
				for(i=0; i<size-1; i++){
					if ( pairedQueryMat->forward_name[skipped][i] == '1' && pairedQueryMat->forward_name[skipped][i-1] == '_'){
						tempname[i] = '2';
					}else{
						tempname[i] = seqname[i];
					}
				}
				tempname[size-1]='\0';
				skipped++;
				if (skipped == numberOfLinesToRead){ break; }
			}*/
			if (skipped == iter){
				memset(pairedQueryMat->reverse_name[iter],'\0',max_readname_length);
				strcpy(pairedQueryMat->reverse_name[iter],seqname);
				memset(seqname,'\0',max_readname_length);
				next=1;
			}else{
				shiftUp(iter,skipped-iter,numberOfLinesToRead);
				next=0;
			}
			/*if ( pairedQueryMat->forward_name[iter+1][0] == '\0' ){
				iter++;
				break;
			}*/
		}else if ( buffer[0] == '@' && whichPair==0){
			for(i=1; i<size; i++){
				if ( buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1]='\0';
			for(i=0; i<max_readname_length; i++){
				singleQueryMat->name[iter][i]='\0';
			}
			strcpy(singleQueryMat->name[iter],seqname);
			for (i=0; i<size-1; i++){
				seqname[i]='\0';
			}
			next=1;
		}else if (next==1){
			for(i=0; i<size; i++){
				query[i]=toupper(buffer[i]);
			}
			query[size]='\0';
			if (whichPair == 0){
				if (opt.reverse_single_read != 1){
					strcpy(singleQueryMat->queryMat[iter],query);
				}else{
					getReverseComplement(query,reverse,max_query_length);
					strcpy(singleQueryMat->queryMat[iter],reverse);
				}
			}
			if (whichPair == 1){
				strcpy(pairedQueryMat->query1Mat[iter],query);
			}
			if (whichPair == 2){
				if (opt.reverse_second_of_paired_read != 1){
					strcpy(pairedQueryMat->query2Mat[iter],query);
				}else{
					getReverseComplement(query,reverse,max_query_length);
					strcpy(pairedQueryMat->query2Mat[iter],reverse);
					for(i=0; i<size; i++){
						reverse[i] = '\0';
					}
				}
			}
			for(i=0; i<size; i++){
				query[i] = '\0';
			}
			iter++;
			next=0;
			if(iter==numberOfLinesToRead){ break; }
		}
	}
	free(buffer);
	free(query);
	free(reverse);
	return iter;
}
//this functions returns the number of 'reads'. This is half the number of actual lines.
//if eof it returns zero.
int readInXNumberOfLines(int numberOfLinesToRead, gzFile query_reads, int whichPair, Options opt, int max_query_length, int max_readname_length){
	char* buffer;
	char* query;
	char* reverse;
	int buffer_size = 0;
	if ( max_query_length > max_readname_length ){
		buffer_size = max_query_length +2;
		buffer = (char *)malloc(sizeof(char)*(max_query_length+2));
	}else{
		buffer_size = max_readname_length +2;
		buffer = (char *)malloc(sizeof(char)*(max_readname_length+2));
	}
	char seqname[max_readname_length];
	int size;
	char *s;
	int i;
	int iter=0;
	int next=0;
	query = (char *)malloc(sizeof(char)*max_query_length+2);
	reverse = (char *)malloc(sizeof(char)*max_query_length+2);
	for(i=0; i<max_query_length+2; i++){
		query[i]='\0';
		reverse[i]='\0';
	}
	int first_iter=1;
	int line_number = 0;
	// Set the current file being processed based on whichPair
	const char* current_filename = NULL;
	if (whichPair == 1) {
		current_filename = opt.read1_file;
	} else if (whichPair == 2) {
		current_filename = opt.read2_file;
	}
	if (current_filename) {
		crash_set_current_file(current_filename);
	}
	
	while(gzgets(query_reads,buffer,buffer_size)!=NULL){
		line_number++;
		crash_set_current_file_line(current_filename, line_number);  // Update current file and line number
		s = strtok(buffer,"\n");
		size = strlen(s);
		if (first_iter==1){
			if ( buffer[0] != '>' ){
				printf("Query reads are not in FASTA format. Try specifying -q if using FASTQ reads.\n");
				exit(-1);
			}
			first_iter=0;
		}
		
		// Check for data corruption patterns
		if (buffer[0] == '>') {
			// Header line - validate format
			if (size < 2) {
				crash_flag_corruption(current_filename, line_number, "Empty FASTA header");
			}
		} else {
			// Sequence line - check for corruption patterns
			if (buffer[0] == ' ' || buffer[0] == '\t') {
				crash_flag_corruption(current_filename, line_number, "Sequence line starts with whitespace");
			}
			if (size > 0 && size < 10) {  // Very short sequences are suspicious
				crash_flag_corruption(current_filename, line_number, "Suspiciously short sequence");
			}
		}
		if ( buffer[0] == '>' && whichPair==1 ){
			for(i=1; i<size; i++){
				if ( buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			memset(pairedQueryMat->forward_name[iter],'\0',max_readname_length);
			strcpy(pairedQueryMat->forward_name[iter],seqname);
			memset(seqname,'\0',max_readname_length);
			next=1;
		}else if (buffer[0] == '>' && whichPair==2) {
			for(i=1; i<size; i++){
				if ( buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			char tempname[max_readname_length];
			memset(tempname,'\0',max_readname_length);
			for(i=0; i<size-1; i++){
				if (pairedQueryMat->forward_name[iter][i] == '1' && pairedQueryMat->forward_name[iter][i-1] == '_'){
					tempname[i] ='2';
				}else{
					tempname[i] = seqname[i];
				}
			}
			tempname[size-1]='\0';
			int skipped=iter;
			/*while(strcmp(tempname,seqname)!=0){
				for(i=0; i<size-1; i++){
					if ( pairedQueryMat->forward_name[skipped][i] == '1' && pairedQueryMat->forward_name[skipped][i-1] == '_'){
						tempname[i] = '2';
					}else{
						tempname[i] = seqname[i];
					}
					tempname[size-1]='\0';
					skipped++;
					if (skipped == numberOfLinesToRead){ break; }
				}
			}*/
			if (skipped == iter){
				memset(pairedQueryMat->reverse_name[iter],'\0',max_readname_length);
				strcpy(pairedQueryMat->reverse_name[iter],seqname);
				memset(seqname,'\0',max_readname_length);
				next=1;
			}else{
				shiftUp(iter,skipped-iter,numberOfLinesToRead);
				next=0;
			}
			/*if ( pairedQueryMat->forward_name[iter+1][0] == '\0' ){
				iter++;
				break;
			}*/
			/*for(i=0; i < MAXREADNAME; i++){
				pairedQueryMat->reverse_name[iter][i]='\0';
			}
			strcpy(pairedQueryMat->reverse_name[iter],seqname);
			for (i=0; i<size-1; i++){
				seqname[i]='\0';
			}*/
		}else if (buffer[0] == '>' && whichPair==0){
			for(i=1; i<size; i++){
				if ( buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			for(i=0; i<max_readname_length; i++){
				singleQueryMat->name[iter][i]='\0';
			}
			strcpy(singleQueryMat->name[iter],seqname);
			for (i=0; i<size-1; i++){
				seqname[i]='\0';
			}
			next=1;
		}else{
			//char query[size];
			//char* query;
			//query = (char *)malloc(sizeof(char)*size);
			for(i=0; i<size; i++){
				query[i]=toupper(buffer[i]);
			}
			query[size]='\0';
			if ( whichPair == 0 ){
				if (opt.reverse_single_read != 1){
					strcpy(singleQueryMat->queryMat[iter],query);
				}else{
					getReverseComplement(query,reverse,max_query_length);
					strcpy(singleQueryMat->queryMat[iter],reverse);
				}
			}
			if ( whichPair == 1 ){
				strcpy(pairedQueryMat->query1Mat[iter],query);
			}
			if ( whichPair == 2 ){
				//char* reverse;
				//reverse = (char *)malloc(sizeof(char)*size);
				if (opt.reverse_second_of_paired_read != 1){
					strcpy(pairedQueryMat->query2Mat[iter],query);
				}else{
					getReverseComplement(query,reverse,max_query_length+2);
					strcpy(pairedQueryMat->query2Mat[iter],reverse);
				//free(reverse);
				//strcpy(pairedQueryMat->query2Mat[iter],query);
					for(i=0; i<size; i++){
						reverse[i] = '\0';
					}
				}
			}
			for(i=0; i<size; i++){
				query[i] = '\0';
			}
			//free(query);
			iter++;
			next=0;
			if(iter==numberOfLinesToRead)
				break;
			}
	}
	free(buffer);
	free(query);
	free(reverse);
	return iter;
}
void shiftUp(int iter, int jump, int numberOfLinesToRead){
	int i;
	for(i=iter; i<numberOfLinesToRead-jump; i++){
		strcpy(pairedQueryMat->forward_name[i],pairedQueryMat->forward_name[i+jump]);
		strcpy(pairedQueryMat->query1Mat[i],pairedQueryMat->query1Mat[i+jump]);
	}
	for(i=numberOfLinesToRead-jump; i<numberOfLinesToRead; i++){
		memset(pairedQueryMat->forward_name[i],'\0',MAXREADNAME);
	}
}
int readReferenceTree( gzFile referenceTree, int* name_specs){
	char buffer[BUFFER_SIZE];
	char acc_name[30];
	int up[2], taxIndex[2];
	char *s;
	int max_lineTaxonomy, max_tax_name, max_nodename, treeNumber, nodeNumber, down, depth, success, i, j, k, firstIter, numberOfTrees;
	firstIter = 1;
	int last_logged_tree = -1;  // Track tree transitions for TREE_LOADED logging
	char* refTreeFlag = buffer;
	while (refTreeFlag != NULL ){
		if ( firstIter == 1 ){
			refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
			if(refTreeFlag == NULL) {
				break;
			}
			s = strtok(buffer, "\n");
			if ( s == NULL ){
				success = 0;
			}else{
				success = sscanf(s, "%d", &numberOfTrees);
			}
			refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
			if(refTreeFlag == NULL) {
				break;
			}
			s = strtok(buffer, "\n");
			if ( s == NULL ){
				success = 0;
			}else{
				success = sscanf(s, "%d", &max_nodename);
				name_specs[0]=max_nodename;
			}
			refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
			s = strtok(buffer, "\n");
			if ( s == NULL ){
				success = 0;
			}else{
				success = sscanf(s, "%d", &max_tax_name);
				name_specs[1]=max_tax_name;
			}
			refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
			s = strtok(buffer, "\n");
			if ( s == NULL ){
				success = 0;
			}else{
				success = sscanf(s, "%d", &max_lineTaxonomy);
				name_specs[2]=max_lineTaxonomy;
			}
			numbaseArr = (int*)malloc(numberOfTrees*(sizeof(int)));
			rootArr = (int*)malloc(numberOfTrees*(sizeof(int)));
			numspecArr= (int*)malloc(numberOfTrees*(sizeof(int)));
			for(i=0; i<numberOfTrees; i++){	
				refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
				if ( refTreeFlag == NULL ){
					break;
				}
				s = strtok(buffer, "\n");
				if ( s == NULL ){
					success = 0;
				}else{
					success = sscanf(s, "%d\t%d\t%d",&(numbaseArr[i]),&(rootArr[i]),&(numspecArr[i]));
				}	
			}
			//for(i=0;i<numberOfTrees;i++){
			//	printf("Tree %d Numbase: %d, Root: %d, Numspec %d\n",i,numbaseArr[i],rootArr[i],numspecArr[i]);
			//}
			allocateMemoryForTaxArr(numberOfTrees,max_tax_name);
			for(i=0;i<numberOfTrees;i++){
				for(j=0; j<numspecArr[i]; j++){
					refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
					if(refTreeFlag == NULL) {
						break;
					}
					s = strtok(buffer,";\n");
					strcpy(taxonomyArr[i][j][0], s ? s : "");
					for(k=1;k<7;k++){
						s = strtok(NULL,";\n");
						strcpy(taxonomyArr[i][j][k], s ? s : "");
					}
				}
			}
			treeArr = malloc(numberOfTrees*sizeof(struct node *));
			for (i=0; i<numberOfTrees; i++){
				allocateTreeArrMemory(i,max_nodename);

				// Per-tree logging at DEBUG level
				LOG_DEBUG("Tree %d allocated: nodes=%d, bases=%d", i, 2*numspecArr[i]-1, numbaseArr[i]);

				// TSV logging if enabled
				if (g_tsv_log_file) {
					resource_stats_t stats;
					get_resource_stats(&stats);
					fprintf(g_tsv_log_file, "%.3f\tTREE_ALLOCATED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d,nodes=%d,bases=%d\n",
							stats.wall_time_sec,
							stats.memory_rss_kb / 1024.0,
							stats.memory_vm_size_kb / 1024.0,
							stats.memory_vm_rss_peak_kb / 1024.0,
							stats.user_time_sec,
							stats.system_time_sec,
							i, 2*numspecArr[i]-1, numbaseArr[i]);
					fflush(g_tsv_log_file);
				}
			}
			//allocatetreememory_for_nucleotide_Arr(numberOfTrees);
			firstIter=0;
		}
		refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
		if (refTreeFlag == NULL){
			break;
		}
		s = strtok(buffer, "\n");
		sscanf(s, "%d %d %d %d %d %d %d %d %s",&treeNumber, &nodeNumber, &(up[0]), &(up[1]), &down, &depth, &(taxIndex[0]), &(taxIndex[1]), acc_name);
		treeArr[treeNumber][nodeNumber].up[0] = up[0];
		treeArr[treeNumber][nodeNumber].up[1] = up[1];
		treeArr[treeNumber][nodeNumber].down = down;
		treeArr[treeNumber][nodeNumber].depth = depth;
		treeArr[treeNumber][nodeNumber].taxIndex[0] = taxIndex[0];
		treeArr[treeNumber][nodeNumber].taxIndex[1] = taxIndex[1];
		if ( up[0] == -1 && up[1] == -1){
			strcpy(treeArr[treeNumber][nodeNumber].name,acc_name);
		}
		for(i=0; i<numbaseArr[treeNumber]; i++){
			refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
			if(refTreeFlag == NULL) {
				break;
			}
			for (j=0; j<4; j++){
				if (j==0){
					s = strtok(buffer,"\t");
				}else{
					s = strtok(NULL,"\t");
				}
				if ( s==NULL ){
					success = 0;
				}else{
					success = sscanf(s, PP_FORMAT, &(treeArr[treeNumber][nodeNumber].posteriornc[PP_IDX(i, j)]));
				}
			}
		}

		// Log when we finish loading a tree's posteriors (tree number changes)
		if (treeNumber != last_logged_tree && last_logged_tree >= 0) {
			LOG_DEBUG("Tree %d posteriors loaded", last_logged_tree);

			if (g_tsv_log_file) {
				resource_stats_t stats;
				get_resource_stats(&stats);
				fprintf(g_tsv_log_file, "%.3f\tTREE_LOADED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d\n",
						stats.wall_time_sec,
						stats.memory_rss_kb / 1024.0,
						stats.memory_vm_size_kb / 1024.0,
						stats.memory_vm_rss_peak_kb / 1024.0,
						stats.user_time_sec,
						stats.system_time_sec,
						last_logged_tree);
				fflush(g_tsv_log_file);
			}
		}
		last_logged_tree = treeNumber;
	}

	// Log the final tree
	if (last_logged_tree >= 0) {
		LOG_DEBUG("Tree %d posteriors loaded (final)", last_logged_tree);

		if (g_tsv_log_file) {
			resource_stats_t stats;
			get_resource_stats(&stats);
			fprintf(g_tsv_log_file, "%.3f\tTREE_LOADED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d\n",
					stats.wall_time_sec,
					stats.memory_rss_kb / 1024.0,
					stats.memory_vm_size_kb / 1024.0,
					stats.memory_vm_rss_peak_kb / 1024.0,
					stats.user_time_sec,
					stats.system_time_sec,
					last_logged_tree);
			fflush(g_tsv_log_file);
		}
	}
	return numberOfTrees;
}

/**
 * Read reference database from binary format (.trkb)
 * Populates the same global structures as readReferenceTree()
 *
 * @param filename Path to .trkb file
 * @param name_specs Output array: [max_nodename, max_tax_name, max_line_taxonomy]
 * @return Number of trees loaded, or -1 on error
 */
int readReferenceBinary(const char *filename, int *name_specs) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        LOG_ERROR("Cannot open binary reference file: %s", filename);
        return -1;
    }

    // === Validate Header ===
    uint8_t magic[4];
    if (fread(magic, 1, 4, fp) != 4) {
        LOG_ERROR("Failed to read magic number");
        fclose(fp);
        return -1;
    }

    if (magic[0] != TRONKO_MAGIC_0 || magic[1] != TRONKO_MAGIC_1 ||
        magic[2] != TRONKO_MAGIC_2 || magic[3] != TRONKO_MAGIC_3) {
        LOG_ERROR("Invalid binary format magic number");
        fclose(fp);
        return -1;
    }

    // Read version
    uint8_t version_major = fgetc(fp);
    uint8_t version_minor = fgetc(fp);
    LOG_DEBUG("Binary format version: %d.%d", version_major, version_minor);

    if (version_major > 1) {
        LOG_ERROR("Unsupported binary format version: %d.%d", version_major, version_minor);
        fclose(fp);
        return -1;
    }

    // Read flags
    uint8_t endianness = fgetc(fp);
    uint8_t precision = fgetc(fp);

    if (endianness != 0x01) {
        LOG_ERROR("Only little-endian binary format is supported");
        fclose(fp);
        return -1;
    }

    if (precision != 0x01) {
        LOG_ERROR("Only float precision binary format is supported");
        fclose(fp);
        return -1;
    }

    // Skip header CRC and reserved
    read_u32(fp);  // header_crc
    read_u32(fp);  // reserved

    // Read section offsets
    uint64_t taxonomy_offset = read_u64(fp);
    uint64_t node_offset = read_u64(fp);
    uint64_t posterior_offset = read_u64(fp);
    uint64_t total_size = read_u64(fp);

    LOG_DEBUG("Section offsets: taxonomy=%lu, nodes=%lu, posteriors=%lu, total=%lu",
              (unsigned long)taxonomy_offset, (unsigned long)node_offset,
              (unsigned long)posterior_offset, (unsigned long)total_size);

    // Skip reserved bytes to reach global metadata
    fseek(fp, BINARY_FILE_HEADER_SIZE, SEEK_SET);

    // === Read Global Metadata (16 bytes) ===
    int32_t numberOfTrees = read_i32(fp);
    int32_t max_nodename = read_i32(fp);
    int32_t max_tax_name = read_i32(fp);
    int32_t max_lineTaxonomy = read_i32(fp);

    name_specs[0] = max_nodename;
    name_specs[1] = max_tax_name;
    name_specs[2] = max_lineTaxonomy;

    LOG_DEBUG("Database: %d trees, max_nodename=%d, max_tax=%d",
              numberOfTrees, max_nodename, max_tax_name);

    // === Allocate Global Arrays ===
    numbaseArr = (int*)malloc(numberOfTrees * sizeof(int));
    rootArr = (int*)malloc(numberOfTrees * sizeof(int));
    numspecArr = (int*)malloc(numberOfTrees * sizeof(int));

    if (!numbaseArr || !rootArr || !numspecArr) {
        LOG_ERROR("Failed to allocate tree metadata arrays");
        fclose(fp);
        return -1;
    }

    // === Read Tree Metadata (12 bytes per tree) ===
    for (int i = 0; i < numberOfTrees; i++) {
        numbaseArr[i] = read_i32(fp);
        rootArr[i] = read_i32(fp);
        numspecArr[i] = read_i32(fp);

        LOG_DEBUG("Tree %d: numbase=%d, root=%d, numspec=%d",
                  i, numbaseArr[i], rootArr[i], numspecArr[i]);
    }

    // === Read Taxonomy Section ===
    fseek(fp, taxonomy_offset, SEEK_SET);

    // Allocate taxonomy array (reuse existing function)
    allocateMemoryForTaxArr(numberOfTrees, max_tax_name);

    for (int t = 0; t < numberOfTrees; t++) {
        uint32_t tree_tax_size = read_u32(fp);
        read_u32(fp);  // reserved
        (void)tree_tax_size;  // Size used for seeking, we read strings directly

        for (int s = 0; s < numspecArr[t]; s++) {
            for (int l = 0; l < 7; l++) {  // 7 taxonomy levels
                uint16_t len = read_u16(fp);
                if (len > 0 && len <= (uint16_t)(max_tax_name + 1)) {
                    if (fread(taxonomyArr[t][s][l], 1, len, fp) != len) {
                        LOG_ERROR("Failed to read taxonomy string");
                        fclose(fp);
                        return -1;
                    }
                } else if (len > 0) {
                    LOG_WARN("Taxonomy string length %d exceeds max %d", len, max_tax_name);
                    // Skip oversized string
                    fseek(fp, len, SEEK_CUR);
                }
            }
        }

        LOG_DEBUG("Tree %d: read %d taxonomy entries", t, numspecArr[t]);
    }

    // === Read Node Section ===
    fseek(fp, node_offset, SEEK_SET);

    // Allocate tree array
    treeArr = malloc(numberOfTrees * sizeof(struct node *));
    if (!treeArr) {
        LOG_ERROR("Failed to allocate tree array");
        fclose(fp);
        return -1;
    }

    for (int t = 0; t < numberOfTrees; t++) {
        uint32_t num_nodes = read_u32(fp);
        int expected_nodes = 2 * numspecArr[t] - 1;

        if ((int)num_nodes != expected_nodes) {
            LOG_ERROR("Node count mismatch for tree %d: got %d, expected %d",
                      t, num_nodes, expected_nodes);
            fclose(fp);
            return -1;
        }

        // Allocate nodes for this tree (reuse existing function pattern)
        allocateTreeArrMemory(t, max_nodename);

        // Store name offsets for later
        uint32_t *name_offsets = calloc(num_nodes, sizeof(uint32_t));
        if (!name_offsets) {
            LOG_ERROR("Failed to allocate name offset array");
            fclose(fp);
            return -1;
        }

        // Read node records (32 bytes each)
        for (int n = 0; n < (int)num_nodes; n++) {
            treeArr[t][n].up[0] = read_i32(fp);
            treeArr[t][n].up[1] = read_i32(fp);
            treeArr[t][n].down = read_i32(fp);
            treeArr[t][n].depth = read_i32(fp);
            treeArr[t][n].taxIndex[0] = read_i32(fp);
            treeArr[t][n].taxIndex[1] = read_i32(fp);
            name_offsets[n] = read_u32(fp);
            read_u32(fp);  // reserved
        }

        // Read name table (for leaf nodes)
        long name_table_start = ftell(fp);
        for (int n = 0; n < (int)num_nodes; n++) {
            if (name_offsets[n] > 0) {
                fseek(fp, name_table_start + name_offsets[n] - 1, SEEK_SET);
                uint16_t len = read_u16(fp);
                if (len > 0 && len <= (uint16_t)(max_nodename + 1)) {
                    // Name already allocated by allocateTreeArrMemory
                    if (fread(treeArr[t][n].name, 1, len, fp) != len) {
                        LOG_ERROR("Failed to read node name");
                        free(name_offsets);
                        fclose(fp);
                        return -1;
                    }
                }
            }
        }

        free(name_offsets);

        LOG_DEBUG("Tree %d: read %d node structures", t, num_nodes);

        // TSV logging if enabled
        if (g_tsv_log_file) {
            resource_stats_t stats;
            get_resource_stats(&stats);
            fprintf(g_tsv_log_file, "%.3f\tTREE_ALLOCATED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d,nodes=%d,bases=%d\n",
                    stats.wall_time_sec,
                    stats.memory_rss_kb / 1024.0,
                    stats.memory_vm_size_kb / 1024.0,
                    stats.memory_vm_rss_peak_kb / 1024.0,
                    stats.user_time_sec,
                    stats.system_time_sec,
                    t, (int)num_nodes, numbaseArr[t]);
            fflush(g_tsv_log_file);
        }
    }

    // === Read Posterior Section ===
    fseek(fp, posterior_offset, SEEK_SET);

    LOG_INFO("Loading posteriors from binary format...");

    for (int t = 0; t < numberOfTrees; t++) {
        int num_nodes = 2 * numspecArr[t] - 1;
        int numbase = numbaseArr[t];

        for (int n = 0; n < num_nodes; n++) {
            // Bulk read all posteriors for this node
            // posteriors array already allocated by allocateTreeArrMemory
            size_t count = numbase * 4;

#ifdef OPTIMIZE_MEMORY
            // Direct read - posteriors are float, format uses float
            if (fread(treeArr[t][n].posteriornc, sizeof(float), count, fp) != count) {
                LOG_ERROR("Failed to read posteriors for tree %d node %d", t, n);
                fclose(fp);
                return -1;
            }
#else
            // Need to convert float to double
            float *temp = malloc(count * sizeof(float));
            if (!temp) {
                LOG_ERROR("Failed to allocate temp buffer for posterior conversion");
                fclose(fp);
                return -1;
            }
            if (fread(temp, sizeof(float), count, fp) != count) {
                LOG_ERROR("Failed to read posteriors for tree %d node %d", t, n);
                free(temp);
                fclose(fp);
                return -1;
            }
            for (size_t i = 0; i < count; i++) {
                treeArr[t][n].posteriornc[i] = (double)temp[i];
            }
            free(temp);
#endif
        }

        LOG_DEBUG("Tree %d: read posteriors for %d nodes", t, num_nodes);

        // TSV logging
        if (g_tsv_log_file) {
            resource_stats_t stats;
            get_resource_stats(&stats);
            fprintf(g_tsv_log_file, "%.3f\tTREE_LOADED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d\n",
                    stats.wall_time_sec,
                    stats.memory_rss_kb / 1024.0,
                    stats.memory_vm_size_kb / 1024.0,
                    stats.memory_vm_rss_peak_kb / 1024.0,
                    stats.user_time_sec,
                    stats.system_time_sec,
                    t);
            fflush(g_tsv_log_file);
        }
    }

    // Skip footer validation (not critical for loading)
    fclose(fp);

    LOG_INFO("Loaded %d trees from binary format", numberOfTrees);

    return numberOfTrees;
}

/**
 * Read reference database from gzipped binary format (.trkb.gz)
 * Populates the same global structures as readReferenceTree()
 *
 * @param filename Path to .trkb.gz file
 * @param name_specs Output array: [max_nodename, max_tax_name, max_line_taxonomy]
 * @return Number of trees loaded, or -1 on error
 */
int readReferenceBinaryGzipped(const char *filename, int *name_specs) {
    gzFile gz = gzopen(filename, "rb");
    if (!gz) {
        LOG_ERROR("Cannot open gzipped binary reference file: %s", filename);
        return -1;
    }

    // === Validate Header ===
    uint8_t magic[4];
    if (gzread(gz, magic, 4) != 4) {
        LOG_ERROR("Failed to read magic number (gzip)");
        gzclose(gz);
        return -1;
    }

    if (magic[0] != TRONKO_MAGIC_0 || magic[1] != TRONKO_MAGIC_1 ||
        magic[2] != TRONKO_MAGIC_2 || magic[3] != TRONKO_MAGIC_3) {
        LOG_ERROR("Invalid binary format magic number (gzip)");
        gzclose(gz);
        return -1;
    }

    // Read version
    uint8_t version_major, version_minor;
    gzread(gz, &version_major, 1);
    gzread(gz, &version_minor, 1);
    LOG_DEBUG("Binary format version: %d.%d (gzip)", version_major, version_minor);

    if (version_major > 1) {
        LOG_ERROR("Unsupported binary format version: %d.%d", version_major, version_minor);
        gzclose(gz);
        return -1;
    }

    // Read flags
    uint8_t endianness, precision;
    gzread(gz, &endianness, 1);
    gzread(gz, &precision, 1);

    if (endianness != 0x01) {
        LOG_ERROR("Only little-endian binary format is supported");
        gzclose(gz);
        return -1;
    }

    if (precision != 0x01) {
        LOG_ERROR("Only float precision binary format is supported");
        gzclose(gz);
        return -1;
    }

    // Skip header CRC and reserved
    gz_read_u32(gz);  // header_crc
    gz_read_u32(gz);  // reserved

    // Read section offsets (we still need these for validation, but we read sequentially)
    uint64_t taxonomy_offset = gz_read_u64(gz);
    uint64_t node_offset = gz_read_u64(gz);
    uint64_t posterior_offset = gz_read_u64(gz);
    uint64_t total_size = gz_read_u64(gz);

    LOG_DEBUG("Section offsets: taxonomy=%lu, nodes=%lu, posteriors=%lu, total=%lu",
              (unsigned long)taxonomy_offset, (unsigned long)node_offset,
              (unsigned long)posterior_offset, (unsigned long)total_size);

    // Skip to global metadata (seek to byte 64)
    gzseek(gz, BINARY_FILE_HEADER_SIZE, SEEK_SET);

    // === Read Global Metadata (16 bytes) ===
    int32_t numberOfTrees = gz_read_i32(gz);
    int32_t max_nodename = gz_read_i32(gz);
    int32_t max_tax_name = gz_read_i32(gz);
    int32_t max_lineTaxonomy = gz_read_i32(gz);

    name_specs[0] = max_nodename;
    name_specs[1] = max_tax_name;
    name_specs[2] = max_lineTaxonomy;

    LOG_DEBUG("Database: %d trees, max_nodename=%d, max_tax=%d (gzip)",
              numberOfTrees, max_nodename, max_tax_name);

    // === Allocate Global Arrays ===
    numbaseArr = (int*)malloc(numberOfTrees * sizeof(int));
    rootArr = (int*)malloc(numberOfTrees * sizeof(int));
    numspecArr = (int*)malloc(numberOfTrees * sizeof(int));

    if (!numbaseArr || !rootArr || !numspecArr) {
        LOG_ERROR("Failed to allocate tree metadata arrays");
        gzclose(gz);
        return -1;
    }

    // === Read Tree Metadata (12 bytes per tree) ===
    for (int i = 0; i < numberOfTrees; i++) {
        numbaseArr[i] = gz_read_i32(gz);
        rootArr[i] = gz_read_i32(gz);
        numspecArr[i] = gz_read_i32(gz);

        LOG_DEBUG("Tree %d: numbase=%d, root=%d, numspec=%d",
                  i, numbaseArr[i], rootArr[i], numspecArr[i]);
    }

    // === Read Taxonomy Section ===
    gzseek(gz, taxonomy_offset, SEEK_SET);

    // Allocate taxonomy array
    allocateMemoryForTaxArr(numberOfTrees, max_tax_name);

    for (int t = 0; t < numberOfTrees; t++) {
        uint32_t tree_tax_size = gz_read_u32(gz);
        gz_read_u32(gz);  // reserved
        (void)tree_tax_size;

        for (int s = 0; s < numspecArr[t]; s++) {
            for (int l = 0; l < 7; l++) {  // 7 taxonomy levels
                uint16_t len = gz_read_u16(gz);
                if (len > 0 && len <= (uint16_t)(max_tax_name + 1)) {
                    if (gzread(gz, taxonomyArr[t][s][l], len) != (int)len) {
                        LOG_ERROR("Failed to read taxonomy string (gzip)");
                        gzclose(gz);
                        return -1;
                    }
                } else if (len > 0) {
                    LOG_WARN("Taxonomy string length %d exceeds max %d", len, max_tax_name);
                    // Skip oversized string by reading and discarding
                    char *skip_buf = malloc(len);
                    if (skip_buf) {
                        gzread(gz, skip_buf, len);
                        free(skip_buf);
                    }
                }
            }
        }

        LOG_DEBUG("Tree %d: read %d taxonomy entries (gzip)", t, numspecArr[t]);
    }

    // === Read Node Section ===
    gzseek(gz, node_offset, SEEK_SET);

    // Allocate tree array
    treeArr = malloc(numberOfTrees * sizeof(struct node *));
    if (!treeArr) {
        LOG_ERROR("Failed to allocate tree array");
        gzclose(gz);
        return -1;
    }

    for (int t = 0; t < numberOfTrees; t++) {
        uint32_t num_nodes = gz_read_u32(gz);
        int expected_nodes = 2 * numspecArr[t] - 1;

        if ((int)num_nodes != expected_nodes) {
            LOG_ERROR("Node count mismatch for tree %d: got %d, expected %d",
                      t, num_nodes, expected_nodes);
            gzclose(gz);
            return -1;
        }

        // Allocate nodes for this tree
        allocateTreeArrMemory(t, max_nodename);

        // Read node records into temporary array first (for name offset handling)
        uint32_t *name_offsets = calloc(num_nodes, sizeof(uint32_t));
        if (!name_offsets) {
            LOG_ERROR("Failed to allocate name offset array");
            gzclose(gz);
            return -1;
        }

        // Read all node records (32 bytes each)
        for (int n = 0; n < (int)num_nodes; n++) {
            treeArr[t][n].up[0] = gz_read_i32(gz);
            treeArr[t][n].up[1] = gz_read_i32(gz);
            treeArr[t][n].down = gz_read_i32(gz);
            treeArr[t][n].depth = gz_read_i32(gz);
            treeArr[t][n].taxIndex[0] = gz_read_i32(gz);
            treeArr[t][n].taxIndex[1] = gz_read_i32(gz);
            name_offsets[n] = gz_read_u32(gz);
            gz_read_u32(gz);  // reserved
        }

        // Read name table - needs special handling for gzip
        // Store current position, then read names sequentially
        z_off_t name_table_start = gztell(gz);

        // For gzip, we can't efficiently seek backwards, so we read names in offset order
        // First, find all non-zero offsets and their indices
        typedef struct { uint32_t offset; int node_idx; } offset_entry_t;
        int name_count = 0;
        for (int n = 0; n < (int)num_nodes; n++) {
            if (name_offsets[n] > 0) name_count++;
        }

        offset_entry_t *sorted_offsets = malloc(name_count * sizeof(offset_entry_t));
        if (!sorted_offsets && name_count > 0) {
            LOG_ERROR("Failed to allocate offset sort array");
            free(name_offsets);
            gzclose(gz);
            return -1;
        }

        int idx = 0;
        for (int n = 0; n < (int)num_nodes; n++) {
            if (name_offsets[n] > 0) {
                sorted_offsets[idx].offset = name_offsets[n];
                sorted_offsets[idx].node_idx = n;
                idx++;
            }
        }

        // Sort by offset (simple bubble sort - name_count is small)
        for (int i = 0; i < name_count - 1; i++) {
            for (int j = 0; j < name_count - i - 1; j++) {
                if (sorted_offsets[j].offset > sorted_offsets[j + 1].offset) {
                    offset_entry_t tmp = sorted_offsets[j];
                    sorted_offsets[j] = sorted_offsets[j + 1];
                    sorted_offsets[j + 1] = tmp;
                }
            }
        }

        // Read names in order
        for (int i = 0; i < name_count; i++) {
            gzseek(gz, name_table_start + sorted_offsets[i].offset - 1, SEEK_SET);
            uint16_t len = gz_read_u16(gz);
            int n = sorted_offsets[i].node_idx;
            if (len > 0 && len <= (uint16_t)(max_nodename + 1)) {
                if (gzread(gz, treeArr[t][n].name, len) != (int)len) {
                    LOG_ERROR("Failed to read node name (gzip)");
                    free(sorted_offsets);
                    free(name_offsets);
                    gzclose(gz);
                    return -1;
                }
            }
        }

        free(sorted_offsets);
        free(name_offsets);

        LOG_DEBUG("Tree %d: read %d node structures (gzip)", t, num_nodes);

        // TSV logging
        if (g_tsv_log_file) {
            resource_stats_t stats;
            get_resource_stats(&stats);
            fprintf(g_tsv_log_file, "%.3f\tTREE_ALLOCATED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d,nodes=%d,bases=%d\n",
                    stats.wall_time_sec,
                    stats.memory_rss_kb / 1024.0,
                    stats.memory_vm_size_kb / 1024.0,
                    stats.memory_vm_rss_peak_kb / 1024.0,
                    stats.user_time_sec,
                    stats.system_time_sec,
                    t, (int)num_nodes, numbaseArr[t]);
            fflush(g_tsv_log_file);
        }
    }

    // === Read Posterior Section ===
    gzseek(gz, posterior_offset, SEEK_SET);

    LOG_INFO("Loading posteriors from gzipped binary format...");

    for (int t = 0; t < numberOfTrees; t++) {
        int num_nodes = 2 * numspecArr[t] - 1;
        int numbase = numbaseArr[t];

        for (int n = 0; n < num_nodes; n++) {
            size_t count = numbase * 4;

#ifdef OPTIMIZE_MEMORY
            // Direct read - posteriors are float
            if (gzread(gz, treeArr[t][n].posteriornc, count * sizeof(float)) != (int)(count * sizeof(float))) {
                LOG_ERROR("Failed to read posteriors for tree %d node %d (gzip)", t, n);
                gzclose(gz);
                return -1;
            }
#else
            // Need to convert float to double
            float *temp = malloc(count * sizeof(float));
            if (!temp) {
                LOG_ERROR("Failed to allocate temp buffer for posterior conversion");
                gzclose(gz);
                return -1;
            }
            if (gzread(gz, temp, count * sizeof(float)) != (int)(count * sizeof(float))) {
                LOG_ERROR("Failed to read posteriors for tree %d node %d (gzip)", t, n);
                free(temp);
                gzclose(gz);
                return -1;
            }
            for (size_t i = 0; i < count; i++) {
                treeArr[t][n].posteriornc[i] = (double)temp[i];
            }
            free(temp);
#endif
        }

        LOG_DEBUG("Tree %d: read posteriors for %d nodes (gzip)", t, num_nodes);

        // TSV logging
        if (g_tsv_log_file) {
            resource_stats_t stats;
            get_resource_stats(&stats);
            fprintf(g_tsv_log_file, "%.3f\tTREE_LOADED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d\n",
                    stats.wall_time_sec,
                    stats.memory_rss_kb / 1024.0,
                    stats.memory_vm_size_kb / 1024.0,
                    stats.memory_vm_rss_peak_kb / 1024.0,
                    stats.user_time_sec,
                    stats.system_time_sec,
                    t);
            fflush(g_tsv_log_file);
        }
    }

    gzclose(gz);

    LOG_INFO("Loaded %d trees from gzipped binary format", numberOfTrees);

    return numberOfTrees;
}

// === Zstd streaming reader implementation ===

typedef struct {
    FILE *file;
    ZSTD_DCtx *dctx;
    uint8_t *in_buf;
    uint8_t *out_buf;
    size_t in_buf_size;
    size_t out_buf_size;
    size_t in_pos;
    size_t in_end;
    size_t out_pos;
    size_t out_end;
    int eof;
} zstd_reader_t;

static zstd_reader_t *zstd_reader_open(const char *filename) {
    zstd_reader_t *r = calloc(1, sizeof(zstd_reader_t));
    if (!r) return NULL;

    r->file = fopen(filename, "rb");
    if (!r->file) { free(r); return NULL; }

    r->dctx = ZSTD_createDCtx();
    if (!r->dctx) { fclose(r->file); free(r); return NULL; }

    r->in_buf_size = ZSTD_DStreamInSize();
    r->out_buf_size = ZSTD_DStreamOutSize();
    r->in_buf = malloc(r->in_buf_size);
    r->out_buf = malloc(r->out_buf_size);

    if (!r->in_buf || !r->out_buf) {
        if (r->in_buf) free(r->in_buf);
        if (r->out_buf) free(r->out_buf);
        ZSTD_freeDCtx(r->dctx);
        fclose(r->file);
        free(r);
        return NULL;
    }

    return r;
}

static int zstd_reader_refill(zstd_reader_t *r) {
    if (r->eof) return 0;

    // Move remaining input data to start
    if (r->in_pos > 0 && r->in_pos < r->in_end) {
        memmove(r->in_buf, r->in_buf + r->in_pos, r->in_end - r->in_pos);
        r->in_end -= r->in_pos;
        r->in_pos = 0;
    } else {
        r->in_pos = 0;
        r->in_end = 0;
    }

    // Read more compressed data
    size_t space = r->in_buf_size - r->in_end;
    if (space > 0) {
        size_t n = fread(r->in_buf + r->in_end, 1, space, r->file);
        r->in_end += n;
        if (n == 0 && r->in_end == 0) {
            r->eof = 1;
            return 0;
        }
    }

    ZSTD_inBuffer in = { r->in_buf, r->in_end, r->in_pos };
    ZSTD_outBuffer out = { r->out_buf, r->out_buf_size, 0 };

    size_t ret = ZSTD_decompressStream(r->dctx, &out, &in);
    if (ZSTD_isError(ret)) {
        LOG_ERROR("ZSTD decompression error: %s", ZSTD_getErrorName(ret));
        return -1;
    }

    r->in_pos = in.pos;
    r->out_pos = 0;
    r->out_end = out.pos;

    return (int)out.pos;
}

static size_t zstd_reader_read(zstd_reader_t *r, void *buf, size_t size) {
    size_t total = 0;
    uint8_t *dst = buf;

    while (total < size) {
        if (r->out_pos >= r->out_end) {
            int ret = zstd_reader_refill(r);
            if (ret <= 0) break;
        }

        size_t avail = r->out_end - r->out_pos;
        size_t need = size - total;
        size_t copy = (avail < need) ? avail : need;

        memcpy(dst + total, r->out_buf + r->out_pos, copy);
        r->out_pos += copy;
        total += copy;
    }

    return total;
}

static void zstd_reader_close(zstd_reader_t *r) {
    if (!r) return;
    ZSTD_freeDCtx(r->dctx);
    free(r->in_buf);
    free(r->out_buf);
    fclose(r->file);
    free(r);
}

// Little-endian read helpers for zstd reader
static uint16_t zstd_read_u16(zstd_reader_t *r) {
    uint8_t b[2];
    if (zstd_reader_read(r, b, 2) != 2) {
        LOG_ERROR("Unexpected end of file reading uint16 (zstd)");
        return 0;
    }
    return (uint16_t)b[0] | ((uint16_t)b[1] << 8);
}

static uint32_t zstd_read_u32(zstd_reader_t *r) {
    uint8_t b[4];
    if (zstd_reader_read(r, b, 4) != 4) {
        LOG_ERROR("Unexpected end of file reading uint32 (zstd)");
        return 0;
    }
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

static int32_t zstd_read_i32(zstd_reader_t *r) {
    return (int32_t)zstd_read_u32(r);
}

static uint64_t zstd_read_u64(zstd_reader_t *r) {
    uint8_t b[8];
    if (zstd_reader_read(r, b, 8) != 8) {
        LOG_ERROR("Unexpected end of file reading uint64 (zstd)");
        return 0;
    }
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) {
        v |= ((uint64_t)b[i] << (i * 8));
    }
    return v;
}

/**
 * Read reference database from zstd-compressed binary format (.trkb)
 * Populates the same global structures as readReferenceTree()
 *
 * @param filename Path to zstd-compressed .trkb file
 * @param name_specs Output array: [max_nodename, max_tax_name, max_line_taxonomy]
 * @return Number of trees loaded, or -1 on error
 */
int readReferenceBinaryZstd(const char *filename, int *name_specs) {
    zstd_reader_t *r = zstd_reader_open(filename);
    if (!r) {
        LOG_ERROR("Cannot open zstd-compressed binary reference file: %s", filename);
        return -1;
    }

    // === Validate Header ===
    uint8_t magic[4];
    if (zstd_reader_read(r, magic, 4) != 4) {
        LOG_ERROR("Failed to read magic number (zstd)");
        zstd_reader_close(r);
        return -1;
    }

    if (magic[0] != TRONKO_MAGIC_0 || magic[1] != TRONKO_MAGIC_1 ||
        magic[2] != TRONKO_MAGIC_2 || magic[3] != TRONKO_MAGIC_3) {
        LOG_ERROR("Invalid binary format magic number (zstd)");
        zstd_reader_close(r);
        return -1;
    }

    // Read version
    uint8_t version_bytes[2];
    zstd_reader_read(r, version_bytes, 2);
    uint8_t version_major = version_bytes[0];
    uint8_t version_minor = version_bytes[1];
    LOG_DEBUG("Binary format version: %d.%d (zstd)", version_major, version_minor);

    if (version_major > 1) {
        LOG_ERROR("Unsupported binary format version: %d.%d", version_major, version_minor);
        zstd_reader_close(r);
        return -1;
    }

    // Read flags
    uint8_t flags[2];
    zstd_reader_read(r, flags, 2);
    uint8_t endianness = flags[0];
    uint8_t precision = flags[1];

    if (endianness != 0x01) {
        LOG_ERROR("Only little-endian binary format is supported");
        zstd_reader_close(r);
        return -1;
    }

    if (precision != 0x01) {
        LOG_ERROR("Only float precision binary format is supported");
        zstd_reader_close(r);
        return -1;
    }

    // Skip header CRC and reserved
    zstd_read_u32(r);  // header_crc
    zstd_read_u32(r);  // reserved

    // Read section offsets (for validation/logging only - we read sequentially)
    uint64_t taxonomy_offset = zstd_read_u64(r);
    uint64_t node_offset = zstd_read_u64(r);
    uint64_t posterior_offset = zstd_read_u64(r);
    uint64_t total_size = zstd_read_u64(r);

    LOG_DEBUG("Section offsets: taxonomy=%lu, nodes=%lu, posteriors=%lu, total=%lu",
              (unsigned long)taxonomy_offset, (unsigned long)node_offset,
              (unsigned long)posterior_offset, (unsigned long)total_size);

    // Skip reserved bytes (16 bytes to complete 64-byte header)
    uint8_t skip[16];
    zstd_reader_read(r, skip, 16);

    // === Read Global Metadata (16 bytes) ===
    int32_t numberOfTrees = zstd_read_i32(r);
    int32_t max_nodename = zstd_read_i32(r);
    int32_t max_tax_name = zstd_read_i32(r);
    int32_t max_lineTaxonomy = zstd_read_i32(r);

    name_specs[0] = max_nodename;
    name_specs[1] = max_tax_name;
    name_specs[2] = max_lineTaxonomy;

    LOG_DEBUG("Database: %d trees, max_nodename=%d, max_tax=%d (zstd)",
              numberOfTrees, max_nodename, max_tax_name);

    // === Allocate Global Arrays ===
    numbaseArr = (int*)malloc(numberOfTrees * sizeof(int));
    rootArr = (int*)malloc(numberOfTrees * sizeof(int));
    numspecArr = (int*)malloc(numberOfTrees * sizeof(int));

    if (!numbaseArr || !rootArr || !numspecArr) {
        LOG_ERROR("Failed to allocate tree metadata arrays");
        zstd_reader_close(r);
        return -1;
    }

    // === Read Tree Metadata (12 bytes per tree) ===
    for (int i = 0; i < numberOfTrees; i++) {
        numbaseArr[i] = zstd_read_i32(r);
        rootArr[i] = zstd_read_i32(r);
        numspecArr[i] = zstd_read_i32(r);

        LOG_DEBUG("Tree %d: numbase=%d, root=%d, numspec=%d",
                  i, numbaseArr[i], rootArr[i], numspecArr[i]);
    }

    // === Read Taxonomy Section ===
    allocateMemoryForTaxArr(numberOfTrees, max_tax_name);

    for (int t = 0; t < numberOfTrees; t++) {
        uint32_t tree_tax_size = zstd_read_u32(r);
        zstd_read_u32(r);  // reserved
        (void)tree_tax_size;

        for (int s = 0; s < numspecArr[t]; s++) {
            for (int l = 0; l < 7; l++) {
                uint16_t len = zstd_read_u16(r);
                if (len > 0 && len <= (uint16_t)(max_tax_name + 1)) {
                    zstd_reader_read(r, taxonomyArr[t][s][l], len);
                } else if (len > 0) {
                    LOG_WARN("Taxonomy string length %d exceeds max %d", len, max_tax_name);
                    // Skip oversized string
                    char *skip_buf = malloc(len);
                    if (skip_buf) {
                        zstd_reader_read(r, skip_buf, len);
                        free(skip_buf);
                    }
                }
            }
        }

        LOG_DEBUG("Tree %d: read %d taxonomy entries (zstd)", t, numspecArr[t]);
    }

    // === Read Node Section ===
    treeArr = malloc(numberOfTrees * sizeof(struct node *));
    if (!treeArr) {
        LOG_ERROR("Failed to allocate tree array");
        zstd_reader_close(r);
        return -1;
    }

    for (int t = 0; t < numberOfTrees; t++) {
        uint32_t num_nodes = zstd_read_u32(r);
        int expected_nodes = 2 * numspecArr[t] - 1;

        if ((int)num_nodes != expected_nodes) {
            LOG_ERROR("Node count mismatch for tree %d: got %d, expected %d",
                      t, num_nodes, expected_nodes);
            zstd_reader_close(r);
            return -1;
        }

        allocateTreeArrMemory(t, max_nodename);

        uint32_t *name_offsets = calloc(num_nodes, sizeof(uint32_t));
        if (!name_offsets) {
            LOG_ERROR("Failed to allocate name offset array");
            zstd_reader_close(r);
            return -1;
        }

        // Read node records
        for (int n = 0; n < (int)num_nodes; n++) {
            treeArr[t][n].up[0] = zstd_read_i32(r);
            treeArr[t][n].up[1] = zstd_read_i32(r);
            treeArr[t][n].down = zstd_read_i32(r);
            treeArr[t][n].depth = zstd_read_i32(r);
            treeArr[t][n].taxIndex[0] = zstd_read_i32(r);
            treeArr[t][n].taxIndex[1] = zstd_read_i32(r);
            name_offsets[n] = zstd_read_u32(r);
            zstd_read_u32(r);  // reserved
        }

        // For zstd, we read names sequentially (can't seek)
        // Names are stored in order after node records
        for (int n = 0; n < (int)num_nodes; n++) {
            if (name_offsets[n] > 0) {
                uint16_t len = zstd_read_u16(r);
                if (len > 0 && len <= (uint16_t)(max_nodename + 1)) {
                    zstd_reader_read(r, treeArr[t][n].name, len);
                }
            }
        }

        free(name_offsets);

        LOG_DEBUG("Tree %d: read %d node structures (zstd)", t, num_nodes);

        if (g_tsv_log_file) {
            resource_stats_t stats;
            get_resource_stats(&stats);
            fprintf(g_tsv_log_file, "%.3f\tTREE_ALLOCATED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d,nodes=%d,bases=%d\n",
                    stats.wall_time_sec,
                    stats.memory_rss_kb / 1024.0,
                    stats.memory_vm_size_kb / 1024.0,
                    stats.memory_vm_rss_peak_kb / 1024.0,
                    stats.user_time_sec,
                    stats.system_time_sec,
                    t, (int)num_nodes, numbaseArr[t]);
            fflush(g_tsv_log_file);
        }
    }

    // === Read Posterior Section ===
    LOG_INFO("Loading posteriors from zstd-compressed binary format...");

    for (int t = 0; t < numberOfTrees; t++) {
        int num_nodes = 2 * numspecArr[t] - 1;
        int numbase = numbaseArr[t];

        for (int n = 0; n < num_nodes; n++) {
            size_t count = numbase * 4;

#ifdef OPTIMIZE_MEMORY
            if (zstd_reader_read(r, treeArr[t][n].posteriornc, count * sizeof(float)) != count * sizeof(float)) {
                LOG_ERROR("Failed to read posteriors for tree %d node %d (zstd)", t, n);
                zstd_reader_close(r);
                return -1;
            }
#else
            float *temp = malloc(count * sizeof(float));
            if (!temp) {
                LOG_ERROR("Failed to allocate temp buffer for posterior conversion");
                zstd_reader_close(r);
                return -1;
            }
            if (zstd_reader_read(r, temp, count * sizeof(float)) != count * sizeof(float)) {
                LOG_ERROR("Failed to read posteriors for tree %d node %d (zstd)", t, n);
                free(temp);
                zstd_reader_close(r);
                return -1;
            }
            for (size_t i = 0; i < count; i++) {
                treeArr[t][n].posteriornc[i] = (double)temp[i];
            }
            free(temp);
#endif
        }

        LOG_DEBUG("Tree %d: read posteriors for %d nodes (zstd)", t, num_nodes);

        if (g_tsv_log_file) {
            resource_stats_t stats;
            get_resource_stats(&stats);
            fprintf(g_tsv_log_file, "%.3f\tTREE_LOADED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d\n",
                    stats.wall_time_sec,
                    stats.memory_rss_kb / 1024.0,
                    stats.memory_vm_size_kb / 1024.0,
                    stats.memory_vm_rss_peak_kb / 1024.0,
                    stats.user_time_sec,
                    stats.system_time_sec,
                    t);
            fflush(g_tsv_log_file);
        }
    }

    zstd_reader_close(r);

    LOG_INFO("Loaded %d trees from zstd-compressed binary format", numberOfTrees);

    return numberOfTrees;
}

int setNumbase_setNumspec(int numberOfPartitions, int* specs){
	int i, maxNumbase, maxNumSpec;
	maxNumbase=0;
	maxNumSpec=0;
	int numspec_total=0;
	for (i=0; i<numberOfPartitions; i++){
		if (maxNumbase < numbaseArr[i]){
			maxNumbase = numbaseArr[i];
		}
	}
	specs[1] = maxNumbase;
	for( i=0; i<numberOfPartitions; i++){
		numspec_total += numspecArr[i];
		if (numspecArr[i] > maxNumSpec){
			maxNumSpec = numspecArr[i];
		}
	}
	specs[0] = maxNumSpec;
	return numspec_total;
}
void find_specs_for_reads(int* specs, gzFile file, int format){
	int max_name_length = specs[0];
	int max_query_length = specs[1];
	char *buffer = (char *)malloc(sizeof(char)*FASTA_MAXLINE);
	char last_name[FASTA_MAXLINE];
	char *s;
	int size = 0;
	while(gzgets(file,buffer,FASTA_MAXLINE) != NULL){
		s = strtok(buffer,"\n");
		if ( s == NULL || s[0] == '\0' ){
			if ( buffer[0] == '>' || buffer[0] == '@' ){
				fprintf(stderr,"Fatal: encountered an empty header at record for read \"%s\" — aborting.\n",last_name[0] ? last_name : "<unknown>");
			}else{
				fprintf(stderr,"Fatal: encountered an empty sequence for read \"%s\" — aborting.\n",last_name[0] ? last_name : "<unknown>");
			}
			free(buffer);
			gzclose(file);
			exit(EXIT_FAILURE);
		}
		if (buffer[0] == '>' || buffer[0] == '@' ){
			strncpy(last_name, s, FASTA_MAXLINE-1);
			last_name[FASTA_MAXLINE-1] = '\0';
			size = strlen(s);
			if (max_name_length < size ){
				max_name_length = size;
			}
		}else{
			size = strlen(s);
			if ( max_query_length < size ){
				max_query_length = size;
			}
		}
	}
	specs[0] = max_name_length;
	specs[1] = max_query_length;
	free(buffer);
}

/**
 * CompressedFile-based version of find_specs_for_reads
 * Supports gzip and zstd compressed FASTA/FASTQ files
 */
void find_specs_for_reads_cf(int* specs, CompressedFile* file, int format){
	int max_name_length = specs[0];
	int max_query_length = specs[1];
	char *buffer = (char *)malloc(sizeof(char)*FASTA_MAXLINE);
	char last_name[FASTA_MAXLINE];
	char *s;
	int size = 0;
	while(cf_gets(buffer, FASTA_MAXLINE, file) != NULL){
		s = strtok(buffer,"\n");
		if ( s == NULL || s[0] == '\0' ){
			if ( buffer[0] == '>' || buffer[0] == '@' ){
				fprintf(stderr,"Fatal: encountered an empty header at record for read \"%s\" — aborting.\n",last_name[0] ? last_name : "<unknown>");
			}else{
				fprintf(stderr,"Fatal: encountered an empty sequence for read \"%s\" — aborting.\n",last_name[0] ? last_name : "<unknown>");
			}
			free(buffer);
			cf_close(file);
			exit(EXIT_FAILURE);
		}
		if (buffer[0] == '>' || buffer[0] == '@' ){
			strncpy(last_name, s, FASTA_MAXLINE-1);
			last_name[FASTA_MAXLINE-1] = '\0';
			size = strlen(s);
			if (max_name_length < size ){
				max_name_length = size;
			}
		}else{
			size = strlen(s);
			if ( max_query_length < size ){
				max_query_length = size;
			}
		}
	}
	specs[0] = max_name_length;
	specs[1] = max_query_length;
	free(buffer);
}

/**
 * CompressedFile-based version of readInXNumberOfLines
 * Reads FASTA format query sequences from gzip or zstd compressed files
 */
int readInXNumberOfLines_cf(int numberOfLinesToRead, CompressedFile* query_reads, int whichPair, Options opt, int max_query_length, int max_readname_length){
	char* buffer;
	char* query;
	char* reverse;
	int buffer_size = 0;
	if ( max_query_length > max_readname_length ){
		buffer_size = max_query_length +2;
		buffer = (char *)malloc(sizeof(char)*(max_query_length+2));
	}else{
		buffer_size = max_readname_length +2;
		buffer = (char *)malloc(sizeof(char)*(max_readname_length+2));
	}
	char seqname[max_readname_length];
	int size;
	char *s;
	int i;
	int iter=0;
	int next=0;
	query = (char *)malloc(sizeof(char)*max_query_length+2);
	reverse = (char *)malloc(sizeof(char)*max_query_length+2);
	for(i=0; i<max_query_length+2; i++){
		query[i]='\0';
		reverse[i]='\0';
	}
	int first_iter=1;
	int line_number = 0;
	// Set the current file being processed based on whichPair
	const char* current_filename = NULL;
	if (whichPair == 1) {
		current_filename = opt.read1_file;
	} else if (whichPair == 2) {
		current_filename = opt.read2_file;
	}
	if (current_filename) {
		crash_set_current_file(current_filename);
	}

	while(cf_gets(buffer, buffer_size, query_reads)!=NULL){
		line_number++;
		crash_set_current_file_line(current_filename, line_number);
		s = strtok(buffer,"\n");
		size = strlen(s);
		if (first_iter==1){
			if ( buffer[0] != '>' ){
				printf("Query reads are not in FASTA format. Try specifying -q if using FASTQ reads.\n");
				exit(-1);
			}
			first_iter=0;
		}

		// Check for data corruption patterns
		if (buffer[0] == '>') {
			if (size < 2) {
				crash_flag_corruption(current_filename, line_number, "Empty FASTA header");
			}
		} else {
			if (buffer[0] == ' ' || buffer[0] == '\t') {
				crash_flag_corruption(current_filename, line_number, "Sequence line starts with whitespace");
			}
			if (size > 0 && size < 10) {
				crash_flag_corruption(current_filename, line_number, "Suspiciously short sequence");
			}
		}
		if ( buffer[0] == '>' && whichPair==1 ){
			for(i=1; i<size; i++){
				if ( buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			memset(pairedQueryMat->forward_name[iter],'\0',max_readname_length);
			strcpy(pairedQueryMat->forward_name[iter],seqname);
			memset(seqname,'\0',max_readname_length);
			next=1;
		}else if (buffer[0] == '>' && whichPair==2) {
			for(i=1; i<size; i++){
				if ( buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			char tempname[max_readname_length];
			memset(tempname,'\0',max_readname_length);
			for(i=0; i<size-1; i++){
				if (pairedQueryMat->forward_name[iter][i] == '1' && pairedQueryMat->forward_name[iter][i-1] == '_'){
					tempname[i] ='2';
				}else{
					tempname[i] = seqname[i];
				}
			}
			tempname[size-1]='\0';
			int skipped=iter;
			if (skipped == iter){
				memset(pairedQueryMat->reverse_name[iter],'\0',max_readname_length);
				strcpy(pairedQueryMat->reverse_name[iter],seqname);
				memset(seqname,'\0',max_readname_length);
				next=1;
			}else{
				shiftUp(iter,skipped-iter,numberOfLinesToRead);
				next=0;
			}
		}else if (buffer[0] == '>' && whichPair==0){
			for(i=1; i<size; i++){
				if ( buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			for(i=0; i<max_readname_length; i++){
				singleQueryMat->name[iter][i]='\0';
			}
			strcpy(singleQueryMat->name[iter],seqname);
			for (i=0; i<size-1; i++){
				seqname[i]='\0';
			}
			next=1;
		}else{
			for(i=0; i<size; i++){
				query[i]=toupper(buffer[i]);
			}
			query[size]='\0';
			if ( whichPair == 0 ){
				if (opt.reverse_single_read != 1){
					strcpy(singleQueryMat->queryMat[iter],query);
				}else{
					getReverseComplement(query,reverse,max_query_length);
					strcpy(singleQueryMat->queryMat[iter],reverse);
				}
			}
			if ( whichPair == 1 ){
				strcpy(pairedQueryMat->query1Mat[iter],query);
			}
			if ( whichPair == 2 ){
				if (opt.reverse_second_of_paired_read != 1){
					strcpy(pairedQueryMat->query2Mat[iter],query);
				}else{
					getReverseComplement(query,reverse,max_query_length+2);
					strcpy(pairedQueryMat->query2Mat[iter],reverse);
					for(i=0; i<size; i++){
						reverse[i] = '\0';
					}
				}
			}
			for(i=0; i<size; i++){
				query[i] = '\0';
			}
			iter++;
			next=0;
			if(iter==numberOfLinesToRead)
				break;
			}
	}
	free(buffer);
	free(query);
	free(reverse);
	return iter;
}

/**
 * CompressedFile-based version of readInXNumberOfLines_fastq
 * Reads FASTQ format query sequences from gzip or zstd compressed files
 */
int readInXNumberOfLines_fastq_cf(int numberOfLinesToRead, CompressedFile* query_reads, int whichPair, Options opt, int max_query_length, int max_readname_length, int first_iter){
	char* buffer;
	char* query;
	char* reverse;
	int buffer_size = 0;
	if ( max_query_length > max_readname_length ){
		buffer_size = max_query_length + 2;
		buffer = (char *)malloc(sizeof(char)*(max_query_length+2));
	}else{
		buffer_size = max_readname_length + 2;
		buffer = (char *)malloc(sizeof(char)*(max_readname_length+2));
	}
	char* s;
	char seqname[max_readname_length+1];
	int size=0;
	int i=0;
	int iter=0;
	int next=0;
	query = (char *)malloc(sizeof(char)*max_query_length+2);
	reverse = (char *)malloc(sizeof(char)*max_query_length+2);
	for(i=0; i<max_query_length+2; i++){
		query[i] = '\0';
		reverse[i] = '\0';
	}
	int line_number_fastq = 0;
	// Set the current file being processed based on whichPair
	const char* current_filename = NULL;
	if (whichPair == 1) {
		current_filename = opt.read1_file;
	} else if (whichPair == 2) {
		current_filename = opt.read2_file;
	}
	if (current_filename) {
		crash_set_current_file(current_filename);
	}

	while(cf_gets(buffer, buffer_size, query_reads)!=NULL){
		line_number_fastq++;
		crash_set_current_file_line(current_filename, line_number_fastq);
		s = strtok(buffer,"\n");
		size = strlen(s);
		if ( first_iter == 1 ){
			if ( buffer[0] != '@' ){
				printf("Query reads are not in FASTQ format\n");
				exit(-1);
			}
			first_iter=0;
		}

		// Check for data corruption patterns in FASTQ
		if (buffer[0] == '@') {
			if (size < 2) {
				crash_flag_corruption(current_filename, line_number_fastq, "Empty FASTQ header");
			}
		} else if (buffer[0] != '+' && buffer[0] != '>') {
			if (buffer[0] == ' ' || buffer[0] == '\t') {
				crash_flag_corruption(current_filename, line_number_fastq, "FASTQ line starts with whitespace");
			}
			if (size > 0 && size < 10) {
				crash_flag_corruption(current_filename, line_number_fastq, "Suspiciously short FASTQ sequence");
			}
		}
		if ( buffer[0] == '@' && whichPair==1){
			for(i=1; i<size; i++){
				if (buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			memset(pairedQueryMat->forward_name[iter],'\0',max_readname_length);
			strcpy(pairedQueryMat->forward_name[iter],seqname);
			memset(seqname,'\0',max_readname_length);
			next=1;
		}else if ( buffer[0] == '@' && whichPair==2){
			for(i=1; i<size; i++){
				if(buffer[i]==' '){ buffer[i]='_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			char tempname[max_readname_length];
			memset(tempname,'\0',max_readname_length);
			for(i=0; i<size-1; i++){
				if ( pairedQueryMat->forward_name[iter][i] == '1' && pairedQueryMat->forward_name[iter][i-1] == '_'){
					tempname[i] = '2';
				}else{
					tempname[i] = seqname[i];
				}
			}
			tempname[size-1]='\0';
			int skipped=iter;
			if (skipped == iter){
				memset(pairedQueryMat->reverse_name[iter],'\0',max_readname_length);
				strcpy(pairedQueryMat->reverse_name[iter],seqname);
				memset(seqname,'\0',max_readname_length);
				next=1;
			}else{
				shiftUp(iter,skipped-iter,numberOfLinesToRead);
				next=0;
			}
		}else if ( buffer[0] == '@' && whichPair==0){
			for(i=1; i<size; i++){
				if ( buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1]='\0';
			for(i=0; i<max_readname_length; i++){
				singleQueryMat->name[iter][i]='\0';
			}
			strcpy(singleQueryMat->name[iter],seqname);
			for (i=0; i<size-1; i++){
				seqname[i]='\0';
			}
			next=1;
		}else if (next==1){
			for(i=0; i<size; i++){
				query[i]=toupper(buffer[i]);
			}
			query[size]='\0';
			if (whichPair == 0){
				if (opt.reverse_single_read != 1){
					strcpy(singleQueryMat->queryMat[iter],query);
				}else{
					getReverseComplement(query,reverse,max_query_length);
					strcpy(singleQueryMat->queryMat[iter],reverse);
				}
			}
			if (whichPair == 1){
				strcpy(pairedQueryMat->query1Mat[iter],query);
			}
			if (whichPair == 2){
				if (opt.reverse_second_of_paired_read != 1){
					strcpy(pairedQueryMat->query2Mat[iter],query);
				}else{
					getReverseComplement(query,reverse,max_query_length);
					strcpy(pairedQueryMat->query2Mat[iter],reverse);
					for(i=0; i<size; i++){
						reverse[i] = '\0';
					}
				}
			}
			for(i=0; i<size; i++){
				query[i] = '\0';
			}
			iter++;
			next=0;
			if(iter==numberOfLinesToRead){ break; }
		}
	}
	free(buffer);
	free(query);
	free(reverse);
	return iter;
}
