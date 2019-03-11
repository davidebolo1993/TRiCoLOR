#include <iostream>
#include <algorithm>
#include <sdsl/suffix_arrays.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <htslib/vcf.h>
#include <htslib/sam.h>


using namespace sdsl;



struct vcfrec
{
	std::string chrName;
	int svStart;
	int svEnd;
	int gen1;
	int gen2;
	
};



char revcomp(char n) {

	switch(n)

	{  

	case 'A':

		return 'T';

	case 'T':

		return 'A';

	case 'G':

		return 'C';

	case 'C':

		return 'G';

	}

	return 0;
}   



void findReadTrim(bam1_t const* rec, int start, int end, int& leftPos, int& rightPos) { // Taken from Alfred

	int rp = rec->core.pos; // reference pointer
	int sp = 0; // sequence pointer
	
	if (start < rp) start = rp; // don't go below the first coordinate
	if (start == end) end = start +1; // if the interval is entire soft-clipped/inserted, add one coordinate to end

	// Parse the CIGAR
	if (start < end) {
			
		uint32_t* cigar = bam_get_cigar(rec);

		for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {

			if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {

				for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {

					++sp;
					++rp;
		  
					if ((leftPos == -1) && (rp >= start)) leftPos = sp;      
				
					if ((rightPos == -1) && (rp >= end)) {
			
						rightPos = sp;
						//return true;
				
					}
				}
			} 

			else {
				
				if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
				  
					rp += bam_cigar_oplen(cigar[i]);
				
				} 

				else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
				  
					sp += bam_cigar_oplen(cigar[i]);
				
				} 

				else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
					
					sp += bam_cigar_oplen(cigar[i]);
				
				} 

				else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {

					//do nothing
				
				} 

				else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
				  
					rp += bam_cigar_oplen(cigar[i]);
				
				} 

				else {
			  
					//return 1;    
				}


				if ((leftPos == -1) && (rp >= start)) leftPos = sp;
		
				if ((rightPos == -1) && (rp >= end)) {
		
					rightPos = sp;
					//return true;
		
				}
			}
		}
	}
  
	else {
	
		//return false;

	}
}



std::string sequencechecker(bam1_t* rec, int start, int end, csa_wt<>& fm_index_ref, std::string seq) {


	bool found=false;
	int app;
	std::string sub;
	
	while(!found) {

		int leftPos=-1;
		int rightPos=-1;
		findReadTrim(rec,start, end, leftPos, rightPos);

		if (seq.substr(leftPos, rightPos - leftPos) == sub) { // if the previous substring is the same of the actual, extending is not useful; sequence is cleared and return an empty string for further evaluation

			sub.clear();
			break;

		}

		else {

			sub = seq.substr(leftPos, rightPos - leftPos);

		}

		app = count(fm_index_ref,sub);

		if (app == 0) {

			found=true;

		}

		else {

			start-=1; // move by 1 bases left. Move just a little, in order to avoid to increase the probability of an error
			end+=1; // move by 1 bases right. Move just a little, in order to avoid to increase the probability of an error
		}
	}

	return sub; 

}



std::string sequencegetter(samFile* samfile1, hts_idx_t* pointer, int head, int start, int end, csa_wt<>& fm_index_ref) {


	hts_itr_t* iter = sam_itr_queryi(pointer, head, start, end);
	bam1_t* rec = bam_init1();
	std::string sequence, seqtr;
	sequence.resize(rec->core.l_qseq);

	while (sam_itr_next(samfile1, iter, rec) >= 0) {

		std::string sequence;
		sequence.resize(rec->core.l_qseq);
		uint8_t* seq = bam_get_seq(rec);

		for (int i = 0; i < rec->core.l_qseq; ++i) {

			if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) { // discard everythig but primary alignments

				continue;

			}

			else {

				sequence[i] = seq_nt16_str[bam_seqi(seq, i)];

				if (rec->core.flag & BAM_FREVERSE) {

					std::transform(sequence.begin(),sequence.end(),sequence.begin(),revcomp); // just in case the sequence is reverse, transform. In principle, all the primary alignments in bam files are forward.

				}

			}
		}

		seqtr=sequencechecker(rec, start, end, fm_index_ref, sequence);
				
	}

	bam_destroy1(rec);
	hts_itr_destroy(iter);


	return seqtr;
}
 


void get_all_fmi(const boost::filesystem::path& root, const std::string& ext, std::vector<boost::filesystem::path>& ret) {


    if(!boost::filesystem::exists(root) || !boost::filesystem::is_directory(root)) return;

    boost::filesystem::recursive_directory_iterator it(root);
    boost::filesystem::recursive_directory_iterator endit;

    while(it != endit) {

        if(boost::filesystem::is_regular_file(*it) && it->path().extension() == ext) ret.push_back(it->path());
        ++it;

    }

}



bool counter(std::string& seq, int& dim, std::vector<csa_wt<>>& fm_index_ill) {

	int i=0;

	bool found=false;
	bool findable=true;


	while (!found && findable) {

		if (i == dim) {

			findable=false;

		}

		else {

			if (count(fm_index_ill[i], seq) >= 1) {

				found=true;

			}


			else {

				i+=1;

			}

		}

	}

	return found;

}



int main(int argc, char* argv[]) {


	boost::posix_time::ptime timer;


	boost::filesystem::path bcfin(argv[1]);
	boost::filesystem::path hap1bam(argv[2]);
	boost::filesystem::path hap2bam(argv[3]);
	boost::filesystem::path fmref(argv[4]); // is a file
	boost::filesystem::path fmill(argv[5]); // is a directory containin illumina .fmi indexes


	std::string ext=".fmi"; // file ending with .fmi. wchih are .fm indexes
	std::vector<boost::filesystem::path> ret; // create a vector to store all the paths to indexes. Difficult 
	get_all_fmi(fmill,ext,ret); // fill ret with paths

	int dim = ret.size();

	if (dim == 0) {

		std::cout << "No fm indexes in folder" << std::endl;
		return 1;

	}

	// print .fm indexes

	std::cout << "FM indexes in folder are:" << std::endl;

	for (int i=0; i<dim; ++i) {

		std::cout << boost::filesystem::path(ret[i]) << std::endl;

	}

	csa_wt<> fm_index_ref; // variable to load ref fmindex
	std::vector<csa_wt<>> fm_index_ill(dim); // vector to load all the .fm indexes. Same size as ret

	timer = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Loading reference FM index" << std::endl;
	load_from_file(fm_index_ref,fmref.string().c_str());
	std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Done" << std::endl;
	
	std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Loading illumina FM indexes" << std::endl;
	
	//fill vector
	
	for(int i=0; i<dim; ++i) {

		load_from_file(fm_index_ill[i], ret[i].string().c_str());
		timer = boost::posix_time::second_clock::local_time();
		std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Loaded " << ret[i] << std::endl;

	}

	std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Done" << std::endl;

	(void)argc;

	
	// preparing to read .vcf file

	vcfrec rec;

	htsFile *bcf = NULL;
	bcf_hdr_t *header = NULL;
	bcf1_t *record = bcf_init();
	bcf = bcf_open(bcfin.string().c_str(), "r");
	header = bcf_hdr_read(bcf);


	int nsvend = 0;
	int* svend = NULL;
	int ngt_arr = 0;
	int *gt = NULL;


	// preparing to read .bam files


	samFile* samfile1 = sam_open(hap1bam.string().c_str(), "r");
	samFile* samfile2 = sam_open(hap2bam.string().c_str(), "r");


	hts_idx_t* idx1 = sam_index_load(samfile1, hap1bam.string().c_str());
	hts_idx_t* idx2 = sam_index_load(samfile2, hap2bam.string().c_str());

	bam_hdr_t* hdr1 = sam_hdr_read(samfile1);
	bam_hdr_t* hdr2 = sam_hdr_read(samfile2);


	int refIndex1,refIndex2;
		
	std::string seq1, seq2;

	int variant=0;
	int valid=0;


	while(bcf_read(bcf, header, record) == 0) {

		rec.chrName = bcf_hdr_id2name(header, record->rid); // get chromosome name
		rec.svStart = record->pos ;  //0 based            
		bcf_get_info_int32(header, record, "END", &svend, &nsvend);
		rec.svEnd = *svend; // 0 based                              
		bcf_get_format_int32(header, record, "GT", &gt, &ngt_arr);
		rec.gen1=bcf_gt_allele(gt[0]);
		rec.gen2=bcf_gt_allele(gt[1]);

		refIndex1 = bam_name2id(hdr1, rec.chrName.c_str());
		refIndex2 = bam_name2id(hdr2, rec.chrName.c_str());

		std::cout << rec.chrName << "\t" << rec.svStart << "\t" << rec.svEnd << "\t" << rec.gen1 << "|" << rec.gen2 << std::endl;


		if ((rec.gen1 == 1) && ((rec.gen2 == -1) || (rec.gen2 == 0))) {

			variant +=1;

			seq1 = sequencegetter(samfile1, idx1, refIndex1, rec.svStart, rec.svEnd, fm_index_ref);

			if (seq1.empty()) {

				continue;
			}

			else {

				std::cout << seq1 << std::endl;

				if (counter(seq1,dim,fm_index_ill)) {

					valid +=1;
					std::cout << "Confirmed in Illumina" << std::endl;

				}

				else {

					std::cout << "Not confirmed in Illumina" << std::endl;

				}


			} 
		}

		else if (((rec.gen1 == -1) || (rec.gen1 == 0)) && ((rec.gen2 == 1) ||(rec.gen2 == 2))) {


			variant +=1;

			seq2 = sequencegetter(samfile2, idx2, refIndex2, rec.svStart, rec.svEnd, fm_index_ref);


			if (seq2.empty()) {

				continue;
			}

			else {

				std::cout << seq2 << std::endl;

				if (counter(seq2,dim,fm_index_ill)) {

					valid +=1;
					std::cout << "Confirmed in Illumina" << std::endl;

				}

				else {

					std::cout << "Not confirmed in Illumina" << std::endl;


				}


			} 

		}


		else if ((rec.gen1 == 1) && (rec.gen2 == 2)) {


			variant +=2;


			seq1 = sequencegetter(samfile1, idx1, refIndex1, rec.svStart, rec.svEnd, fm_index_ref);


			if (seq1.empty()) {

				continue;
			}

			else {

				std::cout << seq1 << std::endl;

				if (counter(seq1,dim,fm_index_ill)) {

					std::cout << "Confirmed in Illumina" << std::endl;

					valid +=1;

				}

				else {

					std::cout << "Not confirmed in Illumina" << std::endl;


				}

			} 


			seq2 = sequencegetter(samfile2, idx2, refIndex2, rec.svStart, rec.svEnd, fm_index_ref);


			if (seq2.empty()) {

				continue;
			}

			else {

				std::cout << seq2 << std::endl;

				if (counter(seq2,dim,fm_index_ill)) {

					std::cout << "Confirmed in Illumina" << std::endl;

					valid +=1;

				}

				else {

					std::cout << "Not confirmed in Illumina" << std::endl;


				}


			} 

		}

		else {

			continue;

		}


	}

	std::cout << variant << std::endl;
	std::cout << valid << std::endl;

	free(svend);
	free(gt);
	bcf_hdr_destroy(header);
	bcf_destroy(record); 
	bcf_close(bcf);
	
	return 0;

} 
