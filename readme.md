# Mining Stable Quasi-Cliques on Temporal Networks


## Environment

Codes run on C++11 or later. When we run codes on linux, -O2 is recommended to optimizate compilation.

## Dataset file format

### Input file format

- The input file should be started with V and E

	-  V: upper limit of node id or higher  

	-  E: the number of edges

- For each event, there are 3 parameters: vertex1, vertex2, time

	- vertex1, vertex2: the vertices of the changing edge

	- time: the time when the event happens

- the time of any (vertex1,vertex2) must be sorted in ascending order

### Output file format

- There are two output files

	- TGRA.txt : store the spent time of TGRA algorithm 

	- res.(dataset) : dataset refers to dataset name, stroe the spent time of B&B algorithm

- No matter how many times BB.cpp is excuted, we can get only one TGRA.txt file. That means TGRA.txt stores all TGRA algorithm result.

- res.(dataset) is different, results of different parameters of one dataset share one dataset. That means one dataset outputs one res.(dataset) file.  

## Compile and run

1. create an excutable application by compiling "BB.cpp", for example:` g++ -o test BB.cpp`

2. run the excutable application with 5 parameters, which are input file name, gamma,  mode, rho, theta, for example:` ./test dblp 0.8 4 0.6 3`

	- input file name: the name of the file describing the tempral graph

	- gamma

	- mode: five pruning algorithms

	- tho

	- theta