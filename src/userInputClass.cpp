//***************************************************************************
//***************************************************************************
/*
 * File: UserInputClass.cpp
 * Description: The functions and classes used to gather user
 * 	input are defined here.
*/
//***************************************************************************
//***************************************************************************
#include "userInputClass.h"

/*
 * Function: readAndIndexBuffer()
 * Purpose: Store contents of an input file as an array of words.
 *
 * Description: This function stores the contents of an input file
 * 	as an array of strings (called words).  A word is either a
 * 	value or the name of the value.  The name of the value will
 * 	directly preced the numerical value of it.
 * 	Comments can be added in a user input file by using the character: #
 * 	at the beginning of the line in the input file.
 */
void UserInput::fileParsingFunctions::readAndIndexBuffer(char *file)
{
	//Read in Buffer*********************************************************
	int size=0;
	ifstream input;
	// Open the file.  ios::ate sets the "cursor" at end of file.
	input.open (file, ios::in | ios::ate);

	if (input.is_open())
	{
		// Retrieve size of file for memory allocation purposes.
		size = input.tellg();
		// Allocate a sufficient amount of memory to the buffer.
		buffer = new char [size+1];
		//! Move the "get" pointer to the beginning of the input file.
		input.seekg (0, ios::beg);
		//! Read in file contents to buffer.
		input.read (buffer, size);	
		//! Close the file.
		input.close();
	}
	//! Display error if the file could not be opened.
	else 
	{
		cout << "FATAL ERROR! Unable to open Input File." << endl;
		exit(1);
	}

	//Add the ending character.
	buffer[size]='\0';

	//******************************************************
	//Index the Buffer**************************************
	int i=0;
	char c=buffer[i];
	numberOfWords=0;

	while (c!='\0') 
	{
		while (c==' ' || c=='\t' || c=='=' || c=='\n') 
		{
			i++;
			c=buffer[i];


		}
		//Remove Comments
		if (c=='#')
		{
			while (c!='\n')
			{
				i++;
				c=buffer[i];
			}
			i++;
			c=buffer[i];
		}
		//Record Word Number
		else if (c!='\0')
		{
			while (c!=' ' && c!='\t' && c!='=' && c!='\n')
			{
				i++;
				c=buffer[i];
			}
			numberOfWords++;
		}
	}

	wordLocation= new int [numberOfWords];
	wordLength= new int [numberOfWords];

	for (int j=0; j<numberOfWords; ++j)
	{
		wordLocation[j]=0;
		wordLength[j]=0;
	}

	i=0;
	int j=0;
	c=buffer[i];
	
	while (c!='\0') 
	{
		while (c==' ' || c=='\t' || c=='=' || c=='\n') 
		{
			i++;
			c=buffer[i];
		}
		//Remove Comments
		if (c=='#')
		{
			while (c!='\n')
			{
				i++;
				c=buffer[i];
			}
			i++;
			c=buffer[i];
		}
		//Record Word
		else if (c!='\0')
		{
			wordLocation[j]=i;
			while (c!=' ' && c!='\t' && c!='=' && c!='\n')
			{
				i++;
				c=buffer[i];
				wordLength[j]++;
			}
			j++;
		}
	}
	//End: Index the Buffer*********************************
	//******************************************************
}

/*
 * Function: getWordFromBuffer()
 * Purpose: retrieve a word from the buffer at a specified index.
 */
void UserInput::fileParsingFunctions::getWordFromBuffer(
	char *word,
	int wordIndex)
{
	int i=0;
	//Store the Word as a character string.
	for (i=0; i<wordLength[wordIndex]; ++i)
	{
		word[i]=buffer[wordLocation[wordIndex]+i];
	}
	word[i]='\0';
}

/*
 * Function: deleteFileParseMemory()
 * Purpose: Free memory used to parse the input file.
 */
void UserInput::fileParsingFunctions::deleteFileParseMemory()
{
	delete [] buffer;
	delete [] wordLocation;
	delete [] wordLength;
	numberOfWords=0;
}

/*
 * Function: getLogFileName()
 * Purpose: Establish the log file name based on the name of the input file.
 *
 * Description: The log file name is stored as a string.
 */
void UserInput::getInputFileNameInfo(char *file,
	string& inputFilePrefix,
	string& logFileName)
{
	int i=0;
	char fileName_tmp[1000];
	while (file[i] != '.')
	{
		//Removes the "" extension from the input file.
		fileName_tmp[i] = file[i];
		i++;
	}
	fileName_tmp[i]='\0';
	inputFilePrefix=fileName_tmp;
	logFileName=inputFilePrefix+".log";
}
