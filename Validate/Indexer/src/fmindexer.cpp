#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <sdsl/suffix_arrays.hpp>

using namespace sdsl;


int main(int argc, char **argv) {

    (void)argc; 


    csa_wt<> fm_index;

    boost::filesystem::path fast(argv[1]);
    boost::filesystem::path fm(argv[2]);
    std::string type(argv[3]);

    if ((type != "fasta") && (type != "fastq")) {

        std::cout << "Accepted type inputs are fasta and fastq" << std::endl;
        return 1;
    }


    boost::filesystem::path modfast(fast.string() + ".tmp");

    boost::posix_time::ptime timer;

    std::ifstream input(fast.string(), std::ios_base::in | std::ios_base::binary);

    if (input.is_open()) {

        char fcode[4];
        input.seekg(0);
        input.read(fcode, 4);

        // first, check if the input is gzipped


        if (((uint8_t)fcode[0] == (uint8_t)0x1f) && ((uint8_t)fcode[1] == (uint8_t)0x8b)) { // check for magic numbers (gzipped input??)

            timer = boost::posix_time::second_clock::local_time();
            std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Preparing ..." << std::endl;

            std::ifstream to_stream(fast.string().c_str(), std::ios_base::in | std::ios_base::binary);

            boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
            inbuf.push(boost::iostreams::gzip_decompressor());
            inbuf.push(to_stream);
            //Convert streambuf to istream
            std::istream instream(&inbuf);
            //Iterate over lines
            std::string line;
            std::ofstream output(modfast.string());
            
            bool first_seq=true;

            if (type == "fasta") { // input is fasta

                while (std::getline(instream, line)) {

                    if(line.empty()) { // exclude blank lines

                        continue;

                    }

                    else if (line[0] == '>') { // this is the line immediately before the sequence

                        continue;

                    }

                    else { // line is not header and is not empty


                        if (!first_seq) { // if not first time that we see a sequence go to new line and get it.
                              
                            output << std::endl << boost::to_upper_copy(line);
                            continue;
                        
                        }

                        else {

                            first_seq=false;
                            output << boost::to_upper_copy(line); // do not go to new line if it is the first time we see the sequence
                            continue;

                        }                        

                    }

                }
               
            }


            else { // input is .fastq

                bool read_next=false;

                while (std::getline(instream, line)) { // same routine


                    if(line.empty()) { // exclude blank lines

                        continue;

                    }

                    else if (line[0] == '@'){ // this is the line immediately before the sequence

                        //std::cout << line << std::endl;
                        read_next=true;
                        continue;

                    }

                    else { // read only if header before

                        if (read_next==true) {

                            if (!first_seq) {
                        
                                output << std::endl << boost::to_upper_copy(line);
                                read_next=false;
                                continue;
                        
                            }

                            else {

                                first_seq=false;
                                read_next=false;
                                output << boost::to_upper_copy(line);
                                continue;

                            }                        

                        }

                    }
               
                }

            }


            output << std::endl;
            output.close();
            to_stream.close();

        }


        else { // input is not gzipped, avoid using decompressors

            timer = boost::posix_time::second_clock::local_time();
            std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Preparing ..." << std::endl;
           
            std::ifstream to_stream(fast.string().c_str(), std::ios_base::in | std::ios_base::binary);
            std::string line;
            std::ofstream output(modfast.string());

            bool first_seq=true;


            if (type == "fasta") { // input is fasta

                while (std::getline(to_stream, line)) {

                    if(line.empty()) { // exclude blank lines

                        continue;

                    }

                    else if (line[0] == '>') { // this is the line immediately before the sequence

                        continue;

                    }

                    else { // line is not header and is not empty


                        if (!first_seq) { // if not first time that we see a sequence go to new line and get it.
                              
                            output << std::endl << boost::to_upper_copy(line);
                            continue;
                        
                        }

                        else {

                            first_seq=false;
                            output << boost::to_upper_copy(line); // do not go to new line if it is the first time we see the sequence
                            continue;

                        }                        

                    }

                }
               
            }


            else { // input is .fastq

                bool read_next=false;

                while (std::getline(to_stream, line)) { // same routine


                    if(line.empty()) { // exclude blank lines

                        continue;

                    }

                    else if (line[0] == '@'){ // this is the line immediately before the sequence

                        //std::cout << line << std::endl;
                        read_next=true;
                        continue;

                    }

                    else { // read only if header before

                        if (read_next==true) {

                            if (!first_seq) {
                        
                                output << std::endl << boost::to_upper_copy(line);
                                read_next=false;
                                continue;
                        
                            }

                            else {

                                first_seq=false;
                                read_next=false;
                                output << boost::to_upper_copy(line);
                                continue;

                            }                        

                        }

                    }
               
                }

            }
            
            output << std::endl;
            output.close();
            to_stream.close();

        }

        input.close();


        timer = boost::posix_time::second_clock::local_time();
        std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "FM-indexing ..." << std::endl;
        construct(fm_index,modfast.string().c_str(), 1); // construct from file, avoid problems with memory ?
        store_to_file(fm_index,fm.string()); // so that it can be re-used
        timer = boost::posix_time::second_clock::local_time();
        std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "FM-index built and saved" << std::endl;

        //boost::filesystem::remove(modfast);


        return 0;

    }



    else {

        return 1; // if any error occur, return non-0

    }

}
