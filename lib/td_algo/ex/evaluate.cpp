/*
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/


#include "output.h"
#include <fstream.h>
#include <post_an.h>

static const char* HELP ="\n\
Evgenii Rudnyi 2000 (c) All rights reserved \n\
http://www.chem.msu.su/~rudnyi \n\
rudnyi@comp.chem.msu.su \n\n\
Using TDLIB to solve the direct problem - \n\
to compute the equilibrium composition and properties \n\n\
    evaluate [-options] modelfile1[.mod] modelfile2[.mod] ...  \n\n\
           options \n\
    -o  base name of the output file (default - modelfile1) \n\
    -h  help \n\
    -l  license \n";

char* LICENSE = "\
Evgenii Rudnyi 2000, (c) All rights reserved \n\n\
Chemistry Department \n\
Moscow State University \n\
119899 Moscow, Russia \n\
http://www.chem.msu.su/~rudnyi \n\
rudnyi@comp.chem.msu.su \n\n\
This program is free software. See files COPYING and FREE_LICENSE. \n\n\
I will be glad if you like this program. Let me know if you find any \n\
bugs. I would also appreciate your comments. \n\n\
You can change the program if you use TDLIB. \n\n\
Disclaimer of warranty: \n\
This program is supplied as is. I disclaim all warranties, \n\
express or implied, including, without limitation, the warranties of \n\
merchantability and of fitness of this program for any purpose. I assume \n\
no liability for damages direct or consequential, which may result from \n\
the use of this program. \n";

int main(int argc, char *argv[])
{
  try
  {
    vec_of vec;
    bool ModelOut = false;
    string basename;
    string model;
    if (argc == 1)
    {
      cout << HELP << endl;
      return 1;
    }
    for (int i = 1; i < argc; i++)
    {
      if (argv[i][0] == '-' || argv[i][0] == '/')
      {
        switch (tolower(argv[i][1]))
        {
          case 'o' :
            basename = argv[i] + 2;
            if (basename.empty() && ++i < argc)
            {
              basename = argv[i];
            }
            if (basename[basename.size() - 1] == '.')
              basename.erase(basename.size() - 1);
            if (basename.empty())
            {
              cout << "null output file" << endl;
              return 1;
            }
            ModelOut = true;
            break;
          case 'h' : case '?' :
            cout << HELP << endl;
            return 0;
          case 'l' :
            cout << LICENSE << endl; 
            return 0;
          default:
            cout << "Option is not known " << endl;
            cout << HELP << endl;
            return 0;
        }
      }
      else
      {
        model = argv[i];
        ifstream in(model.c_str(), ios::in);
        if (!in)
        {
          model += ".mod";
          in.open(model.c_str(), ios::in);
          if (!in)
          {
            cout << "no file - " << model << " - skipped" << endl;
            continue;
          }
        }
        ReadModelFile(in, model, vec);
      }
    }
    if (!ModelOut)
      basename = "con";
    coef::EstimateComputed();
    for (vec_of_ci i = vec.begin(); i != vec.end(); ++i)
      (*i).out(basename);
//    if (ModelOut)
//    {
//      basename += ".mod";
//      ofstream out(basename.c_str());
// does not work, the interface has been changed
//      WriteModelFile(out, vec);
//    }
    return 0;
  }
  catch (gError &t)
  {
    cout << "Error: " << t.message << endl;
  }
  catch (...)
  {
    cout << "Unknown error" << endl;
  }
  return 1;
}
