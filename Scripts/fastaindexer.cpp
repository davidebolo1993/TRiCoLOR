#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

using namespace sdsl;


int main(int argc, char **argv) {

    csa_wt<> fm_index;


    boost::filesystem::path reference(argv[1]);
    boost::filesystem::path modreference(argv[2]);
    boost::filesystem::path fm(argv[3]);


    boost::posix_time::ptime timer;

    std::ifstream input(reference.string());

    if (input.is_open()) {


        timer = boost::posix_time::second_clock::local_time();
        std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Modifying .fasta before indexing ..." << std::endl;

   
        std::string line;
        std::ofstream output;

        output.open(modreference.string());

        bool first_header=true;

        while (std::getline(input, line)) {

            if(line.empty()) { // exclude blank lines

                continue;

            }


            else if (line[0] == '>') {


                if (!first_header) {
                
                    output << std::endl;
                
                }

                else {

                    first_header=false;

                    continue;


                }


            }
       
            else {
           
                output << boost::to_upper_copy(line);
       
            }

        }


        output << '\n';
        output.close();
        input.close();

        timer = boost::posix_time::second_clock::local_time();
        std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "FM-indexing ..." << std::endl;


        construct(fm_index,modreference.string(), 1);
        store_to_file(fm_index,fm.string());


        timer = boost::posix_time::second_clock::local_time();
        std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Done" << std::endl;


        //std::cout << "AT occurs " << count(fm_index,"AT") << " times.\n";


        return 0;

    }


    else {


        return -1;
    }


}
