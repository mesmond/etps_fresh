//***************************************************************************
//***************************************************************************
/*
 * File: userInputClass.h
 * Description: The functions and classes used to gather user
 * 	input are declared here.
*/
//***************************************************************************
//***************************************************************************

#ifndef INCLUDE_USERINPUT_H_
#define INCLUDE_USERINPUT_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;


/*
 * Class: UserInput
 * Purpose: Store and process user input.
 * Description: A UserInput object is able to read in data from an input
 * 	file and search for a desired quantity.
 */
class UserInput
{
	public:
	struct fileParsingFunctions
	{
		char *buffer;
		int numberOfWords;
		int *wordLength;
		int *wordLocation;

		void readAndIndexBuffer(char *file);
		void getWordFromBuffer(char *word, int wordCount);

		void deleteFileParseMemory();

		template <typename T> T readFromFile(char *file, string desiredVariable)
		{
			char word[1000]; //Storage location for parts of buffer
			readAndIndexBuffer(file);

			int wordIndex=0;
			string thisWord;
			string nextWord;

			T value;
			while (wordIndex < numberOfWords-1)
			{
				getWordFromBuffer(word, wordIndex);
				thisWord=word;
				getWordFromBuffer(word, wordIndex+1);
				nextWord=word;

				if (thisWord==desiredVariable)
				{
					stringstream(nextWord) >> value;
					deleteFileParseMemory();
					return value;
				}
			
				wordIndex++;
			} // End Input Parsing While Loop

			cout << "ERROR: Unable to Find the value for \"" << desiredVariable
				<< "\" in input file!"
				<< "\t Defualt Value=" << value << endl;

			deleteFileParseMemory();

			return value;
		}
	} parseInput;

	void getInputFileNameInfo(char *file,
		string& inputFilePrefix,
		string& logFileName);
};

#endif  // INCLUDE_USERINPUT_H_
