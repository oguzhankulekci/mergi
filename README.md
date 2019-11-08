# mergi
Memory Efficient Reference Genome Indexing

Assuming, sdsl-lite is installed on your computer and  'includePATH' and 'libPATH' are the directories where the sdsl-lite headers and libs are maintained, following commmand builds the application mergi.

g++ ./main.cpp -o mergi -O3 -I includePATH -L libPATH -ldivsufsort -lsdsl -ldivsufsort64 -std=c++11
