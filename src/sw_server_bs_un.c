#include "sw_server.h"

//--------------------------------------------------------------------------------------

int apply_sw_bs_un(sw_server_input_t* input, batch_t *batch) {
  LOG_DEBUG("========= APPLY SW BS UNIFIED START =========\n");

  return BS_UN_POST_PAIR_STAGE;

  LOG_DEBUG("------------ APPLY SW BS UN 4NT START -----------------------\n");
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  genome_t *genome1 = input->genome1_p; // genome ACT (G->A)
  genome_t *genome2 = input->genome2_p; // genome ACT (C->T)
  sw_optarg_t *sw_optarg = &input->sw_optarg;

  // fill gaps between seeds
  /*
  {
    char r[1024];
    size_t start = 169312417;
    size_t end = start + 99;
    genome_read_sequence_by_chr_index(r, 0,
    0, &start, &end, genome2);
    printf("+++++++++++++ genome2 = %s \n", r);
    genome_read_sequence_by_chr_index(r, 0,
    0, &start, &end, genome1);
    printf("+++++++++++++ genome1 = %s \n", r);
  }
  */

  LOG_DEBUG("++++++++++ FILL GAPS 0     ++++++++++\n");
  fill_gaps_bs(mapping_batch, sw_optarg, genome1, genome2, 20, 5, 0);
  merge_seed_regions_bs(mapping_batch, 0);
  fill_end_gaps_bs(mapping_batch, sw_optarg, genome1, genome2, 20, 400, 0);
  LOG_DEBUG("++++++++++  END LIST 0     ++++++++++\n");

  LOG_DEBUG("++++++++++ FILL GAPS 1     ++++++++++\n");
  fill_gaps_bs(mapping_batch, sw_optarg, genome2, genome1, 20, 5, 1);
  merge_seed_regions_bs(mapping_batch, 1);
  fill_end_gaps_bs(mapping_batch, sw_optarg, genome2, genome1, 20, 400, 1);
  LOG_DEBUG("++++++++++  END LIST 1     ++++++++++\n");

  // now we can create the alignments
  fastq_read_t *read;
  array_list_t *fq_batch = mapping_batch->fq_batch;
  
  char *match_seq, *match_qual;
  size_t read_index, read_len, match_len, match_start;
  
  cal_t *cal;
  array_list_t *cal_list = NULL;
  size_t num_cals;
  
  seed_region_t *s;
  cigar_code_t *cigar_code;
  cigar_op_t *first_op;

  float score, norm_score, min_score = input->min_score;

  alignment_t *alignment;
  array_list_t *alignment_list;

  char *p, *optional_fields;
  int optional_fields_length, AS;

  array_list_t **mapping_lists;
  size_t num_targets;
  size_t *targets;

  for (int bs_id = 0; bs_id < 2; bs_id++) {

    if (bs_id == 0) {
      mapping_lists = mapping_batch->mapping_lists;
      num_targets = mapping_batch->num_targets;
      targets = mapping_batch->targets;
    } else {
      mapping_lists = mapping_batch->mapping_lists2;
      num_targets = mapping_batch->num_targets2;
      targets = mapping_batch->targets2;
    }

    for (size_t i = 0; i < num_targets; i++) {
      read_index = targets[i];
      read = (fastq_read_t *) array_list_get(read_index, fq_batch);
      
      cal_list = mapping_lists[read_index];
      num_cals = array_list_size(cal_list);
      
      if (num_cals <= 0) continue;
    
      read_len = read->length;
    
      alignment_list = array_list_new(num_cals, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

      // processing each CAL from this read
      for(size_t j = 0; j < num_cals; j++) {

	// get cal and read index
	cal = array_list_get(j, cal_list);
	if (cal->sr_list->size == 0) continue;
	
	s = (seed_region_t *) linked_list_get_first(cal->sr_list);
	cigar_code = (cigar_code_t *) s->info;
	
	norm_score = cigar_code_get_score(read_len, cigar_code);
	score = norm_score * 100; //read_len;
	LOG_DEBUG_F("score = %0.2f\n", norm_score);

	// filter by SW score
	if (norm_score > min_score) {

	  // update cigar and sequence and quality strings
	  cigar_code_update(cigar_code);
	  LOG_DEBUG_F("\tcigar code = %s\n", new_cigar_code_string(cigar_code));
	  match_start = 0;
	  match_len = cigar_code_nt_length(cigar_code); 
	  first_op = cigar_code_get_first_op(cigar_code);
	  match_start = (first_op && first_op->name == 'H' ? first_op->number : 0);
	  
	  match_seq = (char *) malloc((match_len + 1)* sizeof(char));
	  memcpy(match_seq, &read->sequence[match_start], match_len);
	  match_seq[match_len] = 0;
	  
	  match_qual = (char *) malloc((match_len + 1)* sizeof(char));
	  memcpy(match_qual, &read->quality[match_start], match_len);
	  match_qual[match_len] = 0;
	  
	  // set optional fields
	  optional_fields_length = 100;
	  optional_fields = (char *) calloc(optional_fields_length, sizeof(char));
	  
	  p = optional_fields;
	  AS = (int) norm_score * 100;
	
	  sprintf(p, "ASi");
	  p += 3;
	  memcpy(p, &AS, sizeof(int));
	  p += sizeof(int);
	  
	  sprintf(p, "NHi");
	  p += 3;
	  memcpy(p, &num_cals, sizeof(int));
	  p += sizeof(int);
	  
	  sprintf(p, "NMi");
	  p += 3;
	  memcpy(p, &cigar_code->distance, sizeof(int));
	  p += sizeof(int);
	  
	  assert(read->length == cigar_code_nt_length(cigar_code));
	  
	  // create an alignment and insert it into the list
	  alignment = alignment_new();

	  //read_id = malloc(read->length);
	  size_t header_len = strlen(read->id);
	  char *head_id = (char *) malloc(header_len + 1);
	  
	  get_to_first_blank(read->id, header_len, head_id);
	
	  alignment_init_single_end(head_id, match_seq, match_qual, 
				    cal->strand, cal->chromosome_id - 1, cal->start - 1,
				    strdup(new_cigar_code_string(cigar_code)), 
				    cigar_code_get_num_ops(cigar_code), 
				    norm_score * 254, 1, (num_cals > 1),
				    optional_fields_length, optional_fields, alignment);
	  
	  array_list_insert(alignment, alignment_list);

	  LOG_DEBUG_F("creating alignment (bs_id = %i)...\n", bs_id);
	  //alignment_print(alignment);

	}
	// free cigar
	cigar_code_free(cigar_code);
      }
      
      // free the cal list, and update the mapping list with the alignment list
      array_list_free(cal_list, (void *) cal_free);
      mapping_lists[read_index] = alignment_list;
    }
  }

  LOG_DEBUG("========= END OF APPLY SW BS UNIFIED =========\n");
  //return CONSUMER_STAGE;
  return BS_UN_POST_PAIR_STAGE;
}

//--------------------------------------------------------------------------------------
