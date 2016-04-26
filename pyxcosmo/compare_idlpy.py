#!/usr/bin/env python
#
# 2 oct 2014  rene gastaud add the scalar case only for the idl created fits file 
# 14 March 2015 RG add save and missing
# 24 March 2015 RG add exclude
# 29 March 2015 RG add print final status
# 28 April 2015 RG add rrdiff = rdiff where absolute difference is too big. 
#     version=3
# 29 April 2015 RG remove the comput of rdiff with the first element removed
#     version 4
# 15 May   2015 RG if string do not compute differences version 5
# 07 Jan   2016 Version 6, we do not know what's new

import  datetime
import numpy as np

__version__ = 6

def compare_idlpy(funcname, output_idl, output_py, verbose=False, threshold=2e-8, save=True,missing=False, exclude=[]):

    status = True
    mydate = datetime.datetime.now()
    mynames = output_idl.keys()
    print "debug version 0. of compare_idlpy" 
    if(save):
        text_file = open(funcname+mydate.strftime("_%Y_%m_%d")+"_log.txt", "w")
        text_file.write("test "+funcname+" \n")
        text_file.write(mydate.strftime("Date %A %d/%m/%Y \n") )
        text_file.write("threshold  %3.1e"%(threshold)+" \n")
        text_file.write(  "\n")
    i=0
    for name in mynames:
        rdiff_flag=False
        rrdiff_flag=False
        lstatus=True
        i=i+1
        if(save): text_file.write( name+"\n")
        if not (output_py.has_key(name)):
             print name+ " is missing"
             if (missing): status=0
             if(save): text_file.write( name+"  is missing ********* \n \n")
             continue # go to the next iteration
        # 15 May
        if  isinstance(output_idl[name], (str)): # 15 May
            print 'string not tested for ', name, output_idl[name], output_py[name]
            continue # go to the next iteration
        size1 =  np.size(output_idl[name])
        size2 =  np.size(output_py[name])
        #print name, 'size1', size1, 'size2', size2
        if (size1 > 1):
            shape1 =  output_idl[name].shape
            shape2 =  output_py[name].shape
        else:
             shape1 = size1
             shape2 = size2
        
        if (size1 == size2):
            if (size1 > 1):
                reference = (output_idl[name]).reshape(shape2)
            else:
                reference =  output_idl[name]
            diff = output_py[name]-reference
            if (size1 > 1):
                if (reference.min()*reference.max() > 0):
                    rdiff = diff/reference
                    rdiff_flag=True
                # 29 April 2015
                # remove first element
                # better done with rrdiff
                # if (reference[1:].min() > 0):
                #     rdiff = diff[1:]/reference[1:]
                #     rdiff_flag=True
                if not(rdiff_flag):
                    index = np.where(abs(diff) >  threshold)
                    if (index[0].size > 1):
                        # beware division by zero not checked
                        rrdiff = diff[index]/reference[index]
                        rrdiff_flag=True
                        #print 'rrdiff_flag'
                #print 'hoho', name, rdiff_flag
            else:
                if(reference != 0):
                    rdiff = diff/reference
                    rdiff_flag=True

        else:
            print name+" size different ", size1, size2, shape1, shape2
            print ' do nothing'
            #reference =  output_idl[name][0:size2]
            #diff = output_py[name]-reference
            #if (reference.min()*reference.max() > 0):
            #    rdiff = diff/reference
            #    rdiff_flag=True
        #
        #   status 
        if(rrdiff_flag):
            if (size1 > 1):
                if (rrdiff.max() > threshold): lstatus=False
                if (rrdiff.min() < (-threshold)): lstatus=False
            else:
                 if(abs(rrdiff) > threshold): lstatus=False

        if(rdiff_flag) and not(rrdiff_flag):
            if (size1 > 1):
                if (rdiff.max() > threshold): lstatus=False
                if (rdiff.min() < (-threshold)): lstatus=False
            else:
                 if(abs(rdiff) > threshold): lstatus=False

        if not (rdiff_flag) and not(rrdiff_flag):
           if (size1 > 1):
                if (diff.max() > threshold): lstatus=False
                if (diff.min() < (-threshold)): lstatus=False
           else:
               if(abs(diff) > threshold): lstatus=False
                

        if (verbose):
            if (size1 > 1):
                if (rdiff_flag):
                    print  name, 'size1', size1, 'size2', size2,'idl', output_idl[name].min(), output_idl[name].max(),'diff', diff.min(), diff.max(), 'rdiff', rdiff.min(), rdiff.max(), lstatus
                else:
                    print name, 'size1', size1, 'size2', size2,'idl', output_idl[name].min(), output_idl[name].max(),'diff', diff.min(), diff.max(), 'status=', lstatus
                if (rrdiff_flag):
                    print  name, 'r rdiff',  rrdiff.min(), rrdiff.max()

            else:
                if (rdiff_flag):
                    print name, 'size1', size1, 'size2', size2, 'idl',output_idl[name], 'diff', diff, 'rdiff', rdiff, 'status', lstatus
                else:
                    print name, 'size1', size1, 'size2', size2, 'idl', output_idl[name], 'diff', diff, 'status', lstatus
            #print ' '
        #
        if not(name in exclude): status = lstatus & status
        #
        if(save):
            text_file.write( name+" status="+str(lstatus)+" shape idl="+str(shape1)+" shape py="+str(shape2)+"\n")
            if (size1 > 1):
                text_file.write( name+" min max idl \t \t  %3.1e \t  %3.1e"%( (output_idl[name].min()), (output_idl[name].max()) )+"\n")
                text_file.write( name+" min max py  \t \t %3.1e \t  %3.1e"%( (output_py[name].min()), (output_py[name].max()) )+"\n")
                text_file.write( name+" min max (py-idl) \t %3.1e \t  %3.1e"%( diff.min(), diff.max() )+"\n")
                if (rdiff_flag):
                    text_file.write( name+" min max (py-idl)/idl \t %3.1e \t  %3.1e"%( rdiff.min(), rdiff.max() )+"\n")
                if (rrdiff_flag):
                    text_file.write( name+" min max (py-idl)[index]/idl[index] \t %3.1e \t  %3.1e"%( rrdiff.min(), rrdiff.max() )+"\n")
            else:  
                diff = output_idl[name] - output_py[name]
                if (size2 > 1):
                    diff = output_idl[name] - output_py[name][0]
                if (size2 > 1):
                    text_file.write( name+" idl \t %3.1e python \t  %3.1e  (py-idl)  \t  %3.1e"%(output_idl[name], output_py[name][0], diff  )+"\n")
                else:
                    text_file.write( name+" idl  %3.1e \t python %3.1e  \t (py-idl)  %3.1e"%(output_idl[name], output_py[name], diff  )+"\n")
            text_file.write( "\n")
    print "final status="+str(status)
    if(save): 
        text_file.write("final status="+str(status)+"\n")
        text_file.close()

    return status
