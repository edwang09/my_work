//Jianan Wang 800898865
//This is a code to read a txt file that contains a sequence of numbers separated with semicolon, 
//sort it and output it to a txt file named answer.txt at the same directory.


#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
using namespace std;



int main() {
	fstream openfile;
	string filename;
	cout<< "Enter the file name that you want to sort?"<<endl;
	cin>>filename;
	openfile.open(filename);
	vector<int> numbers;	
// read numbers from file to a vector of number.
	if (openfile.is_open()) {
		string line;
		while (getline(openfile, line, ';')) {
			istringstream iss(line);
			int value;
			while (iss>>value) {
				numbers.push_back(value);
			};
		}
	}
//conduct the sorting.
    for (int j=1; j<numbers.size(); j++) {
		int key =numbers[j];
		int i=j-1;
		while (i>=0 && numbers[i] > key) {
			numbers[i+1]=numbers[i];
			i=i-1;
		}
		numbers[i+1]=key;
	}
//output to the answer.txt file
	ofstream outfile;
	outfile.open("answer.txt");
	for(int i=0;i<numbers.size();i++){
		outfile<<numbers[i]<<";";
	}
	outfile.close();
	cout<<" your output is stored in the file named answer.txt"<< endl;
return 0;
}

