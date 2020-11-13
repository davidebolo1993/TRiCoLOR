#include "spoa/spoa.hpp"
#include <iostream>
#include <fstream>


int main(int argc, char** argv) {

    int alignment_type = 1; //no need to change this, as reads from same group will likely be ~ same length and global alignment works fine.

    //the following can be changed by users

    int match_score = atoi(argv[1]); //default to 5
    int mismatch_score = atoi(argv[2]); //default to -4
    int gap_opening_score = atoi(argv[3]); // default to -8
    int gap_extension_score = atoi(argv[4]); // default to -6

    //commented lines are used in https://github.com/davidebolo1993/TREADMILL

    //reference sequence in .fa format

    //std::string seed = argv[5];
    //std::string reference_seq;

    //std::string line;
    //std::ifstream reference_input(seed, std::ios_base::in | std::ios_base::binary);

    //retrive reference sequence from .fa input

    //if (reference_input.is_open()) {

        //while (std::getline(reference_input, line)) {

            //if (line[0] != '>') {
                
                //reference_seq=line;
            
            //}
        
        //}
        
        //reference_input.close();
    
    //}

    //allele sequences in .fa format

    std::string seqs = argv[5];
    std::vector<std::string> sequences = {};
    std::string line;

    std::ifstream allele_input(seqs, std::ios_base::in | std::ios_base::binary);

    if (allele_input.is_open()) {

        while (std::getline(allele_input, line)) {

            if (line[0] != '>') {

                sequences.push_back(line);
            
            }
        
        }
        
        allele_input.close();
    
    }

    int vsize = static_cast<int>(sequences.size());

    auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(alignment_type),(int8_t) match_score, (int8_t) mismatch_score,(int8_t) gap_opening_score, (int8_t) gap_extension_score);
    auto graph = spoa::createGraph();

    // as from Longshot use of SPOA, slightly modified (https://github.com/pjedge/longshot/blob/master/src/poa/poa_func.cpp)
    // seed the POA with a seed sequence (e.g. reference sequence)
    // this is useful e.g. to reduce noise from noisy PacBio sequence reads.
    // we set num_seeds to num_seqs / 2, i.e. the reference sequence has half the "weight" or multiplicity of the noisy reads
    // the first alignment is added to initialize nodes in the graph
    // but each subsequent alignment is going to be the same, and is only to increase the seed sequence's weight.
    // so for each subsequent addition we save computation by just adding the same alignment onto the graph over and over.


     //if (vsize/2 > 0) {
            
        //auto first_seed_aln = alignment_engine->align(reference_seq, graph);
        //graph->add_alignment(first_seed_aln, reference_seq);
        //auto next_seed_aln = alignment_engine->align(reference_seq, graph);
            
        //for (int i = 1; i < vsize/2; i++){
            
            //graph->add_alignment(next_seed_aln, reference_seq);
        
        //}
    
    //}

    // add each of the real sequences (e.g. noisy sequence reads) to the graph

    for (const auto& it: sequences) {
        
        auto alignment = alignment_engine->align(it, graph);
        graph->add_alignment(alignment, it);

    }

    // generate the consensus sequence

    std::string consensus = graph->generate_consensus();
    sequences.clear();

    fprintf(stdout, ">Consensus\n");
    fprintf(stdout, "%s\n", consensus.c_str());

}