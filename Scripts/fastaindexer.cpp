#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <sdsl/suffix_arrays.hpp>



using namespace sdsl;


int main(int argc, char **argv) {

    csa_wt<> fm_index;


    boost::filesystem::path reference(argv[1]);
    boost::filesystem::path modreference(argv[2]);
    boost::filesystem::path fm(argv[3]);


    boost::posix_time::ptime timer;

    std::ifstream input(reference.string(), std::ios_base::in | std::ios_base::binary);

    if (input.is_open()) {

        char fcode[4];
        input.seekg(0);
        input.read(fcode, 4);

    // first, check if the input is gzipped


        if (((uint8_t)fcode[0] == (uint8_t)0x1f) && ((uint8_t)fcode[1] == (uint8_t)0x8b)) { // check for magic numbers

            timer = boost::posix_time::second_clock::local_time();
            std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Preparing for FM indexing, gzipped input ..." << std::endl;

            std::ifstream to_stream(reference.string(), std::ios_base::in | std::ios_base::binary);

            boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
            inbuf.push(boost::iostreams::gzip_decompressor());
            inbuf.push(to_stream);
            //Convert streambuf to istream
            std::istream instream(&inbuf);
            //Iterate over lines
            std::string line;
            std::ofstream output(modreference.string());

            bool first_header=true;

            while (std::getline(instream, line)) {

                if(line.empty()) { // exclude blank lines

                    continue;

                }

                else if (line[0] == '>') { // exclude headers, but just from the second one


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

            output << std::endl;
            output.close();
            to_stream.close();

        }


        else { // input is not gzipped, avoid using decompressors

            timer = boost::posix_time::second_clock::local_time();
            std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Preparing for FM indexing, not gizzipped input ..." << std::endl;
           
            std::ifstream to_stream(reference.string(), std::ios_base::in | std::ios_base::binary);
            std::string line;
            std::ofstream output(modreference.string());

            bool first_header=true;

            while (std::getline(to_stream, line)) {

                if(line.empty()) { // exclude blank lines

                    continue;

                }

                else if (line[0] == '>') {  // exclude headers, but just from the second one


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

            output << std::endl;
            output.close();
            to_stream.close();

        }


        timer = boost::posix_time::second_clock::local_time();
        std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Starting FM-indexing ..." << std::endl;
        construct(fm_index,modreference.string(), 1); // construct from file, avoid problems with memory ?
        store_to_file(fm_index,fm.string()); // so that it can be re-used
        timer = boost::posix_time::second_clock::local_time();
        std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Done" << std::endl;


        return 0;

    }


    else {

        return 1; // if any error occur

    }

}
