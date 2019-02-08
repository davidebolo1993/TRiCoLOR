#include <iostream>
#include <boost/filesystem.hpp>


int main(int argc, char* argv[])
{
  //do something here
  string a = boost::filesystem::canonical('/home/bolognin/').string();
  std::cout << a << std::endl;
  
  return 0;
}
