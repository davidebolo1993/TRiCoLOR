#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>

// boost libraries

#include <boost/filesystem.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

// htslib

#include <htslib/vcf.h>
#include <htslib/sam.h>

// sdsl

#include <sdsl/suffix_arrays.hpp>
using namespace sdsl;



// structure

struct vcfrec

{
	std::string chrName;
	int svStart;
	int svEnd;
	int gen1;
	int gen2;
	
};

struct IO {
	
	boost::filesystem::path vcf;
	boost::filesystem::path bam1;
	boost::filesystem::path bam2;
	boost::filesystem::path fmref;
	boost::filesystem::path fmill;

};



// edit-distance 1 neighbors

bool in_array(const std::string &value, const std::vector<std::string> &array)
{
    
    return std::find(array.begin(), array.end(), value) != array.end();

}


std::vector<std::string> edits(std::string word) {

    std::string myalphabet="ATCG";
    std::vector<std::string> result;
    std::string temp;
           
    for (size_t i = 0; i < word.size(); i++) {
               
        temp = word.substr(0, i) + word.substr(i+1);
       
        if (!in_array(temp, result)) {
           
            result.push_back(temp);
        }
    }


    //for (size_t i = 0; i < word.size() - 1; i++) {
   
        //temp = word.substr(0, i) + word[i+1] + word[i] + word.substr(i+2);
       
        //if (!in_array(temp, result)) {
           
            //result.push_back(temp);
        //}
    //}


    for (size_t i = 0; i < word.size(); i++) {
       
        for (size_t l=0; l < myalphabet.size(); l++) {
                   
            temp = word.substr(0, i) + myalphabet[l] + word.substr(i+1);

            if (!in_array(temp, result)) {
               
                result.push_back(temp);
            }
        }       
    }
           
    for (size_t i = 0; i < word.size() + 1; i++) {
               
        for (size_t l=0; l < myalphabet.size(); l++) {
                   
            temp = word.substr(0, i) + myalphabet[l] + word.substr(i);

            if (!in_array(temp, result)) {
               
                result.push_back(temp);
            }
        }
   
    }
   
    return result;
       
}

// reverse complement

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


// parse cigar

void findReadTrim(bam1_t const* rec, int32_t start, int32_t end, int32_t& leftPos, int32_t& rightPos) {

	int32_t rp = rec->core.pos; // reference pointer
	int32_t sp = 0; // sequence pointer

	if (start == end) {

		end=start+1;

	}

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

	    			}
	    		}
	    	} else {
	
				if (bam_cigar_op(cigar[i]) == BAM_CDEL) {

					rp += bam_cigar_oplen(cigar[i]);

				} else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	  				
	  				sp += bam_cigar_oplen(cigar[i]);

				} else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {

	  				sp += bam_cigar_oplen(cigar[i]);
	
				} else if(bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	  			
	  			// Nothing
				
				} else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	  
	  				rp += bam_cigar_oplen(cigar[i]);
				} 
				
				if ((leftPos == -1) && (rp >= start)) leftPos = sp;
				if ((rightPos == -1) && (rp >= end)) {
	  
	  				rightPos = sp;
				}
      		}
    	}
  	}
}

//get sequence and check occurrence in FM index

std::string sequencechecker(bam1_t* rec, int start, int end, csa_wt<>& fm_index_ref, std::string seq) {


	bool found=false;
	bool findable=true;
	int app;
	std::string sub;
	

	int leftPos=-1;
	int rightPos=-1;
	findReadTrim(rec,start, end, leftPos, rightPos);

	while(!found && findable) {

		sub=seq.substr(leftPos, rightPos - leftPos);
  		
		app = count(fm_index_ref,sub);

		if (app == 0) {

			found=true;

		}

		else {


			if ((leftPos==0) && (rightPos != seq.size())) rightPos +=1;

			else if ((leftPos!= 0) && (rightPos != seq.size())) {

				rightPos +=1;
				leftPos-=1;
			}

			else if ((leftPos != 0) && (rightPos == seq.size())) leftPos -=1;
			
			else {

				findable=false;
				sub.clear();

			}

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

			if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {

				continue;

			}

			else {

				sequence[i] = seq_nt16_str[bam_seqi(seq, i)];

				if (rec->core.flag & BAM_FREVERSE) {

					std::transform(sequence.begin(),sequence.end(),sequence.begin(),revcomp); 

				}

			}
		
			seqtr=sequencechecker(rec, start, end, fm_index_ref, sequence);

		}
				
	}

	bam_destroy1(rec);
	hts_itr_destroy(iter);


	return seqtr;
}
 
// store FM

void get_all_fmi(const boost::filesystem::path& root, const std::string& ext, std::vector<boost::filesystem::path>& ret) {


    if(!boost::filesystem::exists(root) || !boost::filesystem::is_directory(root)) return;

    boost::filesystem::recursive_directory_iterator it(root);
    boost::filesystem::recursive_directory_iterator endit;

    while(it != endit) {

        if(boost::filesystem::is_regular_file(*it) && it->path().extension() == ext) ret.push_back(it->path());
        ++it;

    }

}

// counter

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

// main

int main(int argc, char* argv[]) {


	IO container;

	boost::program_options::options_description generic("Generic options");

	generic.add_options()

		("help,?", "show help message")
		("bamfile1,1", boost::program_options::value<boost::filesystem::path>(&container.bam1), "bamfile for haplotype 1")
		("bamfile2,2", boost::program_options::value<boost::filesystem::path>(&container.bam2), "bamfile for haplotype 2")
		("reference.fmi,ref.fmi", boost::program_options::value<boost::filesystem::path>(&container.fmref), "fm index for reference .fa sequence")
		("illumina.fmi,ill.fmi", boost::program_options::value<boost::filesystem::path>(&container.fmill), "folder containing fm indexes for illumina .fq sequences")
		;
	
	boost::program_options::options_description hidden("Hidden options");

	hidden.add_options()

		("vcf,v", boost::program_options::value<boost::filesystem::path>(&container.vcf), "input vcf/bcf")
		;


	boost::program_options::positional_options_description pos_args;
	pos_args.add("vcf", -1);
	
	boost::program_options::options_description cmdline_options;
	cmdline_options.add(generic).add(hidden);

	boost::program_options::options_description visible_options;
	visible_options.add(generic);

	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
	boost::program_options::notify(vm);

	if ((vm.count("help")) || (!vm.count("vcf"))) {

		std::cout << "Usage: " << argv[0] << " [OPTIONS] input vcf/bcf file" << std::endl;
		std::cout << visible_options << "\n";
		return -1;
	
	}



	boost::posix_time::ptime timer;


	boost::filesystem::path bcfin(container.vcf);
	boost::filesystem::path hap1bam(container.bam1);
	boost::filesystem::path hap2bam(container.bam2);
	boost::filesystem::path fmref(container.fmref); // is a file
	boost::filesystem::path fmill(container.fmill); // is a directory containin one ore more FMI index


	std::string ext=".fmi"; // file ending with .fmi. wchih are .fm indexes
	std::vector<boost::filesystem::path> ret; // create a vector to store all the paths to indexes. 
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
	int toolong=0;

	while(bcf_read(bcf, header, record) == 0) {

		rec.chrName = bcf_hdr_id2name(header, record->rid); // get chromosome name
		rec.svStart = record->pos ;  //0 based            
		bcf_get_info_int32(header, record, "TREND", &svend, &nsvend);
		rec.svEnd = *svend; // 0 based                              
		bcf_get_format_int32(header, record, "GT", &gt, &ngt_arr);
		rec.gen1=bcf_gt_allele(gt[0]);
		rec.gen2=bcf_gt_allele(gt[1]);

		refIndex1 = bam_name2id(hdr1, rec.chrName.c_str());
		refIndex2 = bam_name2id(hdr2, rec.chrName.c_str());

		std::cout << rec.chrName << "\t" << rec.svStart << "\t" << rec.svEnd << "\t" << rec.gen1 << "|" << rec.gen2 << std::endl;

		if ((rec.gen1 == 1) && ((rec.gen2 == -1) || (rec.gen2 == 0))) {

			variant +=1;
			std::cout << "Processing hap1" << std::endl;
			seq1 = sequencegetter(samfile1, idx1, refIndex1, rec.svStart, rec.svEnd, fm_index_ref);

			if (seq1.empty()) {

				continue;
			
			}

			else if (seq1.length() > 124) {

				toolong +=1;
				continue;
			}

			else {

				std::cout << seq1 << std::endl;

				if (counter(seq1,dim,fm_index_ill)) {

					valid +=1;
					std::cout << "Confirmed in Illumina" << std::endl;

				} else {

					std::cout << "Searching in 1-edit neighbors" << std::endl;
					std::vector<std::string>result1 = edits(seq1);
					int innercounter = 0;

    				for (std::vector<std::string>::const_iterator i = result1.begin(); i != result1.end(); ++i) {
   
    					std::string test1(*i);
    					//std::cout << test1 << std::endl;

  						if (counter(test1,dim,fm_index_ill)) {

							valid +=1;
							std::cout << "Confirmed in Illumina" << std::endl;
							innercounter +=1;
							break;

						}
					}

					if (innercounter == 0) {

						std::cout << "Trying also with 2-edit neighbors" << std::endl;
						int innercounter2=0;
						bool found=false;

						for (std::vector<std::string>::const_iterator i = result1.begin(); i != result1.end(); ++i) {

							std::string test1(*i);
							std::vector<std::string>result1_ = edits(test1);

							for (std::vector<std::string>::const_iterator l = result1_.begin(); l != result1_.end(); ++l) {

								std::string test1_(*l);
								//std::cout << test1_ << std::endl;

								if (counter(test1_, dim, fm_index_ill)) {

									valid +=1;
									std::cout << "Confirmed in Illumina" << std::endl;
									innercounter2 +=1;
									found=true;
									break;
								}

							}

							if (found==true) {

								break;
							}
						
						}

						if (innercounter2 == 0) {

							std::cout << "Not confirmed in Illumina" << std::endl;
						}
  					
    				}

				}

			} 
		}

		else if (((rec.gen1 == -1) || (rec.gen1 == 0)) && ((rec.gen2 == 1) ||(rec.gen2 == 2))) {

			std::cout << "Processing hap2" << std::endl;
			variant +=1;
			seq2 = sequencegetter(samfile2, idx2, refIndex2, rec.svStart, rec.svEnd, fm_index_ref);

			if (seq2.empty()) {

				continue;
			
			} 

			else if (seq2.length() > 124) {

				toolong +=1;
				continue;
			}

			else {

				std::cout << seq2 << std::endl;

				if (counter(seq2,dim,fm_index_ill)) {

					valid +=1;
					std::cout << "Confirmed in Illumina" << std::endl;

				}


				else {

					std::cout << "Searching in 1-edit neighbors" << std::endl;
					std::vector<std::string>result2 = edits(seq2);
					int innercounter = 0;

    				for (std::vector<std::string>::const_iterator i = result2.begin(); i != result2.end(); ++i) {
   
    					std::string test2(*i);
						//std::cout << test2 << std::endl;

  						if (counter(test2,dim,fm_index_ill)) {

							valid +=1;
							std::cout << "Confirmed in Illumina" << std::endl;
							innercounter +=1;
							break;

						}
					}

					if (innercounter == 0) {

						std::cout << "Trying also with 2-edit neighbors" << std::endl;
						int innercounter2=0;
						bool found=false;

						for (std::vector<std::string>::const_iterator i = result2.begin(); i != result2.end(); ++i) {

							std::string test2(*i);
							std::vector<std::string>result2_ = edits(test2);

							for (std::vector<std::string>::const_iterator l = result2_.begin(); l != result2_.end(); ++l) {

								std::string test2_(*l);
								//std::cout << test2_ << std::endl;

								if (counter(test2_, dim, fm_index_ill)) {

									valid +=1;
									std::cout << "Confirmed in Illumina" << std::endl;
									innercounter2 +=1;
									found=true;
									break;
								}

							}

							if (found==true) {

								break;
							}

						}

						if (innercounter2 == 0) {

							std::cout << "Not confirmed in Illumina" << std::endl;
						}
  					
    				}

				}

			} 

		}


		else if ((rec.gen1 == 1) && (rec.gen2 == 2)) {


			variant +=2;
			std::cout << "Processing hap1 and hap2" << std::endl;
			seq1 = sequencegetter(samfile1, idx1, refIndex1, rec.svStart, rec.svEnd, fm_index_ref);

			if (seq1.empty()) {

				continue;
			}


			else if (seq1.length() > 124) {

				toolong +=1;
				continue;
			}

			else {

				std::cout << seq1 << std::endl;

				if (counter(seq1,dim,fm_index_ill)) {

					std::cout << "Confirmed in Illumina" << std::endl;
					valid +=1;

				} else {

					std::cout << "Searching in 1-edit neighbors" << std::endl;
					std::vector<std::string>result1 = edits(seq1);
					int innercounter = 0;

    				for (std::vector<std::string>::const_iterator i = result1.begin(); i != result1.end(); ++i) {
   
    					std::string test1(*i);
    					//std::cout << test1 << std::endl;

  						if (counter(test1,dim,fm_index_ill)) {

							valid +=1;
							std::cout << "Confirmed in Illumina" << std::endl;
							innercounter +=1;
							break;
						}
					}

					if (innercounter == 0) {

						std::cout << "Trying also with 2-edit neighbors" << std::endl;
						int innercounter2=0;
						bool found=false;

						for (std::vector<std::string>::const_iterator i = result1.begin(); i != result1.end(); ++i) {

							std::string test1(*i);
							std::vector<std::string>result1_ = edits(test1);

							for (std::vector<std::string>::const_iterator l = result1_.begin(); l != result1_.end(); ++l) {

								std::string test1_(*l);
								//std::cout << test1_ << std::endl;

								if (counter(test1_, dim, fm_index_ill)) {

									valid +=1;
									std::cout << "Confirmed in Illumina" << std::endl;
									innercounter2 +=1;
									found=true;
									break;
								}

							}

							if (found==true) {

								break;
							}

						
						}

						if (innercounter2 == 0) {

							std::cout << "Not confirmed in Illumina" << std::endl;
						}
  					
    				}

				}

			} 

			seq2 = sequencegetter(samfile2, idx2, refIndex2, rec.svStart, rec.svEnd, fm_index_ref);

			if (seq2.empty()) {

				continue;
			
			}

			else if (seq2.length() > 124) {

				toolong +=1;
				continue;
			}

			else {

				std::cout << seq2 << std::endl;

				if (counter(seq2,dim,fm_index_ill)) {

					std::cout << "Confirmed in Illumina" << std::endl;
					valid +=1;

				}

				else {

					std::cout << "Searching in 1-edit neighbors" << std::endl;
					std::vector<std::string>result2 = edits(seq2);
					int innercounter = 0;

    				for (std::vector<std::string>::const_iterator i = result2.begin(); i != result2.end(); ++i) {
   
    					std::string test2(*i);
    					//std::cout<< test2 << std::endl;

  						if (counter(test2,dim,fm_index_ill)) {

  							//std::cout << test2 << std::endl;
							valid +=1;
							std::cout << "Confirmed in Illumina" << std::endl;
							innercounter +=1;
							break;

						}
					}

					if (innercounter == 0) {

						std::cout << "Trying also with 2-edit neighbors" << std::endl;
						int innercounter2=0;
						bool found=false;

						for (std::vector<std::string>::const_iterator i = result2.begin(); i != result2.end(); ++i) {

							std::string test2(*i);
							std::vector<std::string>result2_ = edits(test2);

							for (std::vector<std::string>::const_iterator l = result2_.begin(); l != result2_.end(); ++l) {

								std::string test2_(*l);
								//std::cout << test2_ << std::endl;

								if (counter(test2_, dim, fm_index_ill)) {

									valid +=1;
									std::cout << "Confirmed in Illumina" << std::endl;
									innercounter2 +=1;
									found=true;
									break;
								}

							}

							if (found==true) {

								break;
							}
						
						}

						if (innercounter2 == 0) {

							std::cout << "Not confirmed in Illumina" << std::endl;
						}
  					
    				}

				}

			} 

		}

		else {

			continue;

		}

	}


	std::ofstream myfile;
	myfile.open ("results.txt");
	myfile << "Total variants: " << variant << std::endl;
	myfile << "Validated variants: " << valid << std::endl;
	myfile << "Variants too long for Illumina: " << toolong << std::endl;
	myfile.close();


	free(svend);
	free(gt);
	bcf_hdr_destroy(header);
	bcf_destroy(record); 
	bcf_close(bcf);
	
	return 0;

}
