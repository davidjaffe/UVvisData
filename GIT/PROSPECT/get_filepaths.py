#!/usr/bin/env python
'''
lone class 
20170103
'''
import sys,os

class get_filepaths():
    def get_filepaths(self,directory,exclude=None):
        '''
        20160906 taken from http://stackoverflow.com/questions/3207219/how-to-list-all-files-of-a-directory-in-python
        This function will generate the file names in a directory 
        tree by walking the tree either top-down or bottom-up. For each 
        directory in the tree rooted at directory top (including top itself), 
        it yields a 3-tuple (dirpath, dirnames, filenames).

        20170103 Modified to allow exclusion based on substring
        '''
        file_paths = []  # List which will store all of the full filepaths.

        # Walk the tree.
        for root, directories, files in os.walk(directory):
            for filename in files:
                # Join the two strings in order to form the full filepath.
                filepath = os.path.join(root, filename)
                Accept = True
                if exclude is not None:
                    if exclude in filename : Accept = False
                if Accept : file_paths.append(filepath)  # Add it to the list.

        return file_paths  # Self-explanatory.
