#include "spoa/spoa.hpp"
#include <iostream>
#include <fstream>


int main(int argc, char** argv) {

    std::vector<std::string> sequences = {};
    std::vector<std::string> headers = {};
    std::string infile = argv[6];

    std::ifstream input(infile, std::ios_base::in | std::ios_base::binary);
    std::string line;
    bool first_seq=true;


    if (input.is_open()) {
        

        while (std::getline(input, line)) {

            if (line[0] == '>') {
                headers.push_back(line);
                continue;
            }

            else {


                if (!first_seq) {
                    
                    sequences.push_back(line);
                }

            else {

                    first_seq=false;
                    sequences.push_back(line);
                }                        

            }

        }

        input.close();

    }


    for (std::vector<std::string>::const_iterator i = headers.begin(); i != headers.end(); ++i) std::cerr << *i << '\n';


    auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(atoi(argv[1])),
        atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    auto graph = spoa::createGraph();

    for (const auto& it: sequences) {
        auto alignment = alignment_engine->align(it, graph);
        graph->add_alignment(alignment, it);
    }

    std::string consensus = graph->generate_consensus();

    fprintf(stdout, ">Consensus\n");
    fprintf(stdout, "%s\n", consensus.c_str());

}
