#include <iostream>
#include <fstream>

// boost libraries

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

// sdsl

#include <sdsl/suffix_arrays.hpp>
using namespace sdsl;

// structure

struct IO {
	
	boost::filesystem::path fasta;
	boost::filesystem::path outfile;
	std::string type;
};

// main

int main(int argc, char **argv) {

	IO container;

	// command line options

	boost::program_options::options_description generic("Generic options");

	generic.add_options()

		("help,?", "show help message")
		("output,o", boost::program_options::value<boost::filesystem::path>(&container.outfile)->default_value("file.fmi"), "output file")
		("type,t", boost::program_options::value<std::string>(&container.type), "type of sequence [fasta,fastq]")
		;
	
	boost::program_options::options_description hidden("Hidden options");

	hidden.add_options()

		("input,i", boost::program_options::value<boost::filesystem::path>(&container.fasta), "input file")
		;

	boost::program_options::positional_options_description pos_args;
	pos_args.add("input", -1);	
	boost::program_options::options_description cmdline_options;
	cmdline_options.add(generic).add(hidden);
	boost::program_options::options_description visible_options;
	visible_options.add(generic);
	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
	boost::program_options::notify(vm);
	
	if ((vm.count("help")) || (!vm.count("input"))) {

		std::cout << "Usage: " << argv[0] << " [OPTIONS] input file" << std::endl;
		std::cout << visible_options << "\n";
		return -1;
	
	}

	// code

	csa_wt<> fm_index;

	boost::filesystem::path fast(container.fasta);
	boost::filesystem::path fm(container.outfile);
	std::string type(container.type);

	if ((type != "fasta") && (type != "fastq")) {

		std::cout << "Accepted type inputs are fasta and fastq" << std::endl;
		return -1;
	
	}

	boost::filesystem::path modfast(fm.string() + ".tmp");
	boost::posix_time::ptime timer;

	if ( !boost::filesystem::exists(modfast.string())) {

		std::ifstream input(fast.string(), std::ios_base::in | std::ios_base::binary);

		if (input.is_open()) {

			char fcode[4];
			input.seekg(0);
			input.read(fcode, 4);

			if (((uint8_t)fcode[0] == (uint8_t)0x1f) && ((uint8_t)fcode[1] == (uint8_t)0x8b)) {

				timer = boost::posix_time::second_clock::local_time();
				std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Preparing ..." << std::endl;
				std::ifstream to_stream(fast.string().c_str(), std::ios_base::in | std::ios_base::binary);
				boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
				inbuf.push(boost::iostreams::gzip_decompressor());
				inbuf.push(to_stream);
				std::istream instream(&inbuf);
				std::string line;
				std::ofstream output(modfast.string());				
				bool first_seq=true;

				if (type == "fasta") {

					while (std::getline(to_stream, line)) {

						if(line.empty()) {

							continue;

						}

						else if (line[0] == '>') {

							if (first_seq) {

								first_seq=false;
								continue;

							}

							else {

								output << std::endl;

							}

						}

						else {


							output << boost::to_upper_copy(line);
						}

					}
				   
				}

				else {

					bool first_seq=true;
					bool first_line=true;
					bool second_line=false;
					bool third_line = false;
					bool fourth_line= false;

					while (std::getline(instream, line)) {

						if(line.empty()) {

							continue;

						}

						else {

							if (first_line==true) {

								first_line=false;
								second_line=true;
								continue;
							
							}

							else if (second_line == true) {

								second_line=false;
								third_line=true;

								if (first_seq == true) {

									first_seq=false;
									output << boost::to_upper_copy(line);

								}

								else {

									output << std::endl << boost::to_upper_copy(line);

								}

								continue;
							
							}

							else if (third_line == true) {

								third_line = false;
								fourth_line=true;
								continue;
							
							}

							else if (fourth_line == true) {

								fourth_line=false;
								first_line=true;
								continue;

							}

						}

					}

				}

				output << std::endl;
				output.close();
				to_stream.close();

			}

			else {

				timer = boost::posix_time::second_clock::local_time();
				std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "Preparing ..." << std::endl;		   
				std::ifstream to_stream(fast.string().c_str(), std::ios_base::in | std::ios_base::binary);
				std::string line;
				std::ofstream output(modfast.string());
				bool first_seq=true;

				if (type == "fasta") {

					while (std::getline(to_stream, line)) {

						if(line.empty()) {

							continue;

						}

						else if (line[0] == '>') {

							if (first_seq) {

								first_seq=false;
								continue;

							}

							else {

								output << std::endl;

							}

						}

						else {


							output << boost::to_upper_copy(line);
						}

					}
				   
				}

				else {

					bool first_seq=true;
					bool first_line=true;
					bool second_line=false;
					bool third_line = false;
					bool fourth_line= false;

					while (std::getline(to_stream, line)) {

						if(line.empty()) {

							continue;

						}

						else {

							if (first_line==true) {

								first_line=false;
								second_line=true;
								continue;
							
							}

							else if (second_line == true) {

								second_line=false;
								third_line=true;

								if (first_seq == true) {

									first_seq=false;
									output << boost::to_upper_copy(line);

								}

								else {

									output << std::endl << boost::to_upper_copy(line);

								}

								continue;
							
							}

							else if (third_line == true) {

								third_line = false;
								fourth_line=true;
								continue;
							
							}

							else if (fourth_line == true) {

								fourth_line=false;
								first_line=true;
								continue;

							}

						}

					}

				}
				
				output << std::endl;
				output.close();
				to_stream.close();

			}

			input.close();
		
		}

		else {

			std::cout  << "Cannot open " << fast.string() << std::endl;
			return -1;
		}
	}

	timer = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "fm-indexing ..." << std::endl;
	construct(fm_index,modfast.string().c_str(), 1);
	store_to_file(fm_index,fm.string());
	timer = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(timer) << "] " << "fm-index built and saved" << std::endl;
	//boost::filesystem::remove(modfast);

	return 0;

}
