#!/usr/bin/env python

# -*- coding: utf8 -*-

'''
cortado __version__ = "1.0" 
https://github.com/staciawyman/cortado
cortado is a reimplementation of CRISPResso https://github.com/lucapinello/CRISPResso
Reimplemented by Stacia Wyman from CRISPResso  Luca Pinello 2015 version 1.0.8
Software pipeline for the analysis of CRISPR-Cas9 genome editing outcomes from deep sequencing data
'''



import sys
import errno
import os
import subprocess as sb
import argparse
import re
import gzip
from collections import defaultdict
import multiprocessing as mp
import pickle as cp
#import cPickle as cp
import unicodedata


import logging
logging.basicConfig(level=logging.INFO,
                     format='%(levelname)-5s @ %(asctime)s:\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info


####Support functions###

_ROOT = os.path.abspath(os.path.dirname(__file__))

def get_data(path):
        return os.path.join(_ROOT, 'data', path)

def check_library(library_name):
        try:
                return __import__(library_name)
        except:
                error('You need to install %s module to use cortado.' % library_name)
                sys.exit(1)

def output_error(msg):
         if os.path.isfile('outputsummary.txt'):
             fh_summary = open('outputsummary.txt',"a")
             append = True
         else:
             fh_summary = open('outputsummary.txt',"w")
             append = False
             fh_summary.write('Sample\tTotalReads\tMergedReads\tPercentMerged\tAlignedReads\tPercentAligned\tUnmodified\t%Unmodified\tCutsiteSubs\tNon-cutsiteSubs\t%CutsiteSubs\tNHEJ\t%NHEJ\n')
         fh_summary.write("%s: ERROR: %s\n" % (args.name,msg))
         fh_summary.close()


def which(program):
        import os
        def is_exe(fpath):
                return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
                if is_exe(program):
                        return program
        else:
                for path in os.environ["PATH"].split(os.pathsep):
                        path = path.strip('"')
                        exe_file = os.path.join(path, program)
                        if is_exe(exe_file):
                                return exe_file

        return None

def check_program(binary_name,download_url=None):
        if not which(binary_name):
                error('You need to install and have the command #####%s##### in your PATH variable to use cortado!\n Please read the documentation!' % binary_name)
                if download_url:
                        error('You can download it from here:%s' % download_url)
                sys.exit(1)


def check_file(filename):
    try:
        with open(filename): pass
    except IOError:
        raise Exception('I cannot open the file: '+filename)

def force_symlink(src, dst):

    if os.path.exists(dst) and os.path.samefile(src,dst):
        return

    try:
        os.symlink(src, dst)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            os.remove(dst)
            os.symlink(src, dst)

nt_complement=dict({'A':'T','C':'G','G':'C','T':'A','N':'N','_':'_','-':'-'})

def reverse_complement(seq):
        return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])

def find_wrong_nt(sequence):
    return list(set(sequence.upper()).difference(set(['A','T','C','G','N'])))


def get_ids_reads_to_remove(fastq_filename,min_bp_quality=20,min_single_bp_quality=0):
    ids_to_remove=set()
    if fastq_filename.endswith('.gz'):
        fastq_handle=gzip.open(fastq_filename)
    else:
        fastq_handle=open(fastq_filename)

    for record in SeqIO.parse(fastq_handle, "fastq"):
        if np.array(record.letter_annotations["phred_quality"]).mean()<min_bp_quality \
        or np.array(record.letter_annotations["phred_quality"]).min()<min_single_bp_quality:
            ids_to_remove.add(record.id)

    return ids_to_remove


def filter_pe_fastq_by_qual(fastq_r1,fastq_r2,output_filename_r1=None,output_filename_r2=None,min_bp_quality=20,min_single_bp_quality=0):

    ids_to_remove_s1=get_ids_reads_to_remove(fastq_r1,min_bp_quality=min_bp_quality,min_single_bp_quality=min_single_bp_quality)
    ids_to_remove_s2=get_ids_reads_to_remove(fastq_r2,min_bp_quality=min_bp_quality,min_single_bp_quality=min_single_bp_quality)

    ids_to_remove=ids_to_remove_s1.union(ids_to_remove_s2)

    if fastq_r1.endswith('.gz'):
        fastq_handle_r1=gzip.open(fastq_r1)
    else:
        fastq_handle_r1=open(fastq_r1)

    if fastq_r2.endswith('.gz'):
        fastq_handle_r2=gzip.open(fastq_r2)
    else:
        fastq_handle_r2=open(fastq_r2)

    if not output_filename_r1:
        output_filename_r1=fastq_r1.replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz'

    if not output_filename_r2:
        output_filename_r2=fastq_r2.replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz'

    #we cannot use with on gzip with python 2.6 :(
    try:
        fastq_filtered_outfile_r1=gzip.open(output_filename_r1,'w+')

        for record in SeqIO.parse(fastq_handle_r1, "fastq"):
            if not record.id in ids_to_remove:
                fastq_filtered_outfile_r1.write(record.format('fastq'))
    except:
        raise Exception('Error handling the fastq_filtered_outfile_r1')

    try:
        fastq_filtered_outfile_r2=gzip.open(output_filename_r2,'w+')

        for record in SeqIO.parse(fastq_handle_r2, "fastq"):
            if not record.id in ids_to_remove:
                fastq_filtered_outfile_r2.write(record.format('fastq'))
    except:
        raise Exception('Error handling the fastq_filtered_outfile_r2')


    return output_filename_r1,output_filename_r2


def filter_se_fastq_by_qual(fastq_filename,output_filename=None,min_bp_quality=20,min_single_bp_quality=0):

        if fastq_filename.endswith('.gz'):
                fastq_handle=gzip.open(fastq_filename)
        else:
                fastq_handle=open(fastq_filename)

        if not output_filename:
                output_filename=fastq_filename.replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz'

        try:
            fastq_filtered_outfile=gzip.open(output_filename,'w+')

            for record in SeqIO.parse(fastq_handle, "fastq"):
                if np.array(record.letter_annotations["phred_quality"]).mean()>=min_bp_quality \
                and np.array(record.letter_annotations["phred_quality"]).min()>=min_single_bp_quality:
                    fastq_filtered_outfile.write(record.format('fastq'))
        except:
                raise Exception('Error handling the fastq_filtered_outfile')


        return output_filename

def get_avg_read_lenght_fastq(fastq_filename):
     cmd=('z' if fastq_filename.endswith('.gz') else '' ) +('cat < %s' % fastq_filename)+\
                  r''' | awk 'BN {n=0;s=0;} NR%4 == 2 {s+=length($0);n++;} END { printf("%d\n",s/n)}' '''
     p = sb.Popen(cmd, shell=True,stdout=sb.PIPE)
     return int(p.communicate()[0].strip())

def get_n_reads_fastq(fastq_filename):
     p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < %s | wc -l" % fastq_filename , shell=True,stdout=sb.PIPE)
     return int(float(p.communicate()[0])/4.0)

matplotlib=check_library('matplotlib')
from matplotlib import font_manager as fm
font = {'size'   : 22}
matplotlib.rc('font', **font)
matplotlib.use('Agg')

plt=check_library('pylab')
from matplotlib import font_manager as fm
from matplotlib import colors as colors_mpl
import matplotlib.gridspec as gridspec

pd=check_library('pandas')
np=check_library('numpy')
Bio=check_library('Bio')

check_program('java')
check_program('flash')
check_program('needle')

sns=check_library('seaborn')
sns.set_context('poster')
sns.set(font_scale=2.2)
sns.set_style('white')

from Bio import SeqIO,pairwise2
#########################################



###EXCEPTIONS############################
class ReorientException(Exception):
    pass

class FlashException(Exception):
    pass

class TrimmomaticException(Exception):
    pass

class NeedleException(Exception):
    pass

class NoReadsAlignedException(Exception):
    pass

class DonorSequenceException(Exception):
    pass

class AmpliconEqualDonorException(Exception):
    pass

class CoreDonorSequenceNotContainedException(Exception):
    pass

class CoreDonorSequenceNotUniqueException(Exception):
    pass

class SgRNASequenceException(Exception):
    pass

class NTException(Exception):
    pass

class ExonSequenceException(Exception):
    pass

class DuplicateSequenceIdException(Exception):
    pass

class NoReadsAfterQualityFiltering(Exception):
    pass

#########################################
def classify_read (row,substitution_positions,insertion_positions_flat,deletion_positions_flat):
    global global_sub_count
    global partial_correction
    global seq_error
    #global ref_donor_diffs
    NHEJ=False
    HDR=False
    PARTIAL=False
    CUTSUB=False
    MIXED=False
    SEQ_ERROR=False

    ref_donor_diffs=[i for i in xrange(len(args.expected_hdr_amplicon_seq)) if args.expected_hdr_amplicon_seq[i] != args.amplicon_seq[i]]
    if include_idxs.intersection(insertion_positions_flat) or  include_idxs.intersection(deletion_positions_flat):
        NHEJ = True

    for i in range(len(substitution_positions)):
        if substitution_positions[i] not in ref_donor_diffs:
            if substitution_positions[i] in include_idxs:
                CUTSUB = True
                #info('===>CUTSUB')
                global_sub_count += 1
            else:
                seq_error += 1
                SEQ_ERROR = True
                #info('===>SEQ_ERROR')
                           

    last = len(ref_donor_diffs)
    count_edits = 0
    for i in range(len(ref_donor_diffs)):
        if ref_donor_diffs[i] in substitution_positions \
        and args.expected_hdr_amplicon_seq[ref_donor_diffs[i]] == row.align_seq[ref_donor_diffs[i]+row.read_offset]:
            if args.main_site == i:
                HDR = True
            if not NHEJ:
                hdr_vector[i]+=1
                count_edits+=1

    # if the number of hdr sites is the total expected, increment counter
    if count_edits == last:
        hdr_vector[last]+=1
    else:
        if not NHEJ and not HDR and count_edits > 0:
            partial_correction+=1
            PARTIAL = True
    
    if NHEJ and HDR: # Main site edit with indel at cutsite
        fh_mixed.write("%s\n%s\n%s\n\n" % (row.ref_seq,row.align_str,row.align_seq))
        fh_indel.write("%s\n%s\n%s\n\n" % (row.ref_seq,row.align_str,row.align_seq))
        return 'NHEJ'
    elif NHEJ: # Indel at cutsite
        fh_indel.write("%s\n%s\n%s\n\n" % (row.ref_seq,row.align_str,row.align_seq))
        return 'NHEJ'
    elif HDR and not CUTSUB: # Main site corrected, no NHEJ, but maybe SUBS? AHA, that's why NHEJ if sub at hdr site...
        fh_hdr.write("%s\n%s\n%s\n\n" % (row.ref_seq,row.align_str,row.align_seq))
        return 'HDR'
    elif CUTSUB:
        fh_cutsub.write("%s\n%s\n%s\n\n" % (row.ref_seq,row.align_str,row.align_seq))
        return 'CUTSUB'
    elif PARTIAL:
        fh_partial.write("%s\n%s\n%s\n\n" % (row.ref_seq,row.align_str,row.align_seq))
        return 'PARTIAL'
    else:
        fh_ref.write("%s\n%s\n%s\n\n" % (row.ref_seq,row.align_str,row.align_seq))
        return 'UNMODIFIED'

def process_df_chunk(df_needle_alignment_chunk):

     MODIFIED_FRAMESHIFT=0
     MODIFIED_NON_FRAMESHIFT=0
     NON_INDELS_NON_FRAMESHIFT=0
     SPLICING_SITES_MODIFIED=0

     #INITIALIZATIONS
     if args.coding_seq:
         PERFORM_FRAMESHIFT_ANALYSIS=True
     else:
         PERFORM_FRAMESHIFT_ANALYSIS=False

     effect_vector_insertion=np.zeros(len_amplicon)
     effect_vector_deletion=np.zeros(len_amplicon)
     effect_vector_mutation=np.zeros(len_amplicon)
     effect_vector_any=np.zeros(len_amplicon)

     effect_vector_insertion_mixed=np.zeros(len_amplicon)
     effect_vector_deletion_mixed=np.zeros(len_amplicon)
     effect_vector_mutation_mixed=np.zeros(len_amplicon)

     effect_vector_insertion_hdr=np.zeros(len_amplicon)
     effect_vector_deletion_hdr=np.zeros(len_amplicon)
     effect_vector_mutation_hdr=np.zeros(len_amplicon)

     effect_vector_insertion_noncoding=np.zeros(len_amplicon)
     effect_vector_deletion_noncoding=np.zeros(len_amplicon)
     effect_vector_mutation_noncoding=np.zeros(len_amplicon)

     hist_inframe=defaultdict(lambda :0)
     hist_frameshift=defaultdict(lambda :0)

     avg_vector_del_all=np.zeros(len_amplicon)
     avg_vector_ins_all=np.zeros(len_amplicon)

     # one or more -+ and \.+
     re_find_indels=re.compile("(-*-)")
     re_find_substitutions=re.compile("(\.*\.)")


     for idx_row,row in df_needle_alignment_chunk.iterrows():

                 #GET THE MUTATIONS POSITIONS
                 # THIS CAN'T EVER BE TRUE??
                 if row.UNMODIFIED:
                     fh_ref.write("%s\n%s\n%s\n\n" % (row.ref_seq,row.align_str,row.align_seq))
                     continue

                 if PERFORM_FRAMESHIFT_ANALYSIS:
                    lenght_modified_positions_exons=[]
                    current_read_exons_modified=False
                    current_read_spliced_modified=False

                 #quantify substitution
                 #GET SUBSTITUTION POSITIONS FROM READ
                 substitution_positions=[]
                 if not args.ignore_substitutions:
                     for p in re_find_substitutions.finditer(row.align_str):
                         st,en=p.span()
                         if args.expected_hdr_amplicon_seq:
                             #if row.aln_donor_seq[st:en] == row.align_seq[st:en]:
                             substitution_positions.append(row.ref_positions[st:en])
                         else: 
                             substitution_positions.append(row.ref_positions[st:en])
                     if substitution_positions:
                         substitution_positions=list(np.hstack(substitution_positions))

                 #GET DELETION POSITIONS FROM READ
                 deletion_positions=[]
                 deletion_positions_flat=[]
                 deletion_sizes=[]
                 if not args.ignore_deletions:
                     for p in re_find_indels.finditer(row.align_seq):
                         st,en=p.span()
                         deletion_positions.append(row.ref_positions[st:en])
                         deletion_sizes.append(en-st)
                     if deletion_positions:
                         deletion_positions_flat=np.hstack(deletion_positions)

                 #quantify insertion
                 insertion_positions=[]
                 insertion_sizes=[]
                 insertion_positions_flat=[]
                 if not args.ignore_insertions:
                     for p in re_find_indels.finditer(row.ref_seq):
                         st,en=p.span()
                         #ref_st=row.ref_positions[st-1] # we report the base preceding the insertion
                         #insertion_positions.append(ref_st)
                         insertion_positions.append([row['ref_positions'][max(0,st-1)],row['ref_positions'][min(len(row['ref_positions'])-1,en)]])
                         insertion_sizes.append(en-st)

                     if insertion_positions:
                         insertion_positions_flat=np.hstack(insertion_positions)
                         #print "%s adding to insertions: %s include: %s" % (insertion_sizes,insertion_positions_flat,include_idxs)



                 ########CLASSIFY READ
                 #hdr  WE HAVE THE DONOR SEQUENCE
                 if args.expected_hdr_amplicon_seq:
                    read_type=classify_read(row,substitution_positions,insertion_positions_flat,deletion_positions_flat)
                    df_needle_alignment_chunk.ix[idx_row,read_type]=True
                 #NON hdr NO DONOR SEQUENCE PROVIDED
                 else:
                    #NHEJ
                    if include_idxs.intersection(insertion_positions_flat) or \
                        include_idxs.intersection(deletion_positions_flat):
                        df_needle_alignment_chunk.ix[idx_row,'NHEJ']=True
                        fh_indel.write("%s\n%s\n%s\n\n" % (row.ref_seq,row.align_str,row.align_seq))
                    if include_idxs.intersection(substitution_positions): 
                        df_needle_alignment_chunk.ix[idx_row,'CUTSUB']=True
                        fh_cutsub.write("%s\n%s\n%s\n\n" % (row.ref_seq,row.align_str,row.align_seq))

                    #UNMODIFIED Can have substitution and be unmodified
                    if not df_needle_alignment_chunk.ix[idx_row,'CUTSUB'] and not df_needle_alignment_chunk.ix[idx_row,'NHEJ']: 
                        df_needle_alignment_chunk.ix[idx_row,'UNMODIFIED']=True
                        fh_ref.write("%s\n%s\n%s\n\n" % (row.ref_seq,row.align_str,row.align_seq))

                 ###CREATE AVERAGE SIGNALS, HERE WE SHOW EVERYTHING...
                 if df_needle_alignment_chunk.ix[idx_row,'MIXED']:
                    effect_vector_mutation_mixed[substitution_positions]+=1
                    effect_vector_deletion_mixed[deletion_positions_flat]+=1
                    effect_vector_insertion_mixed[insertion_positions_flat]+=1

                 elif df_needle_alignment_chunk.ix[idx_row,'HDR']:
                    effect_vector_mutation_hdr[substitution_positions]+=1
                    effect_vector_deletion_hdr[deletion_positions_flat]+=1
                    effect_vector_insertion_hdr[insertion_positions_flat]+=1

                 elif df_needle_alignment_chunk.ix[idx_row,'NHEJ'] and not args.hide_mutations_outside_window_NHEJ:
                    effect_vector_mutation[substitution_positions]+=1
                    effect_vector_deletion[deletion_positions_flat]+=1
                    effect_vector_insertion[insertion_positions_flat]+=1

                 any_positions=np.unique(np.hstack([deletion_positions_flat,insertion_positions_flat,substitution_positions])).astype(int)
                 effect_vector_any[any_positions]+=1


     fh_hdr.close()
     fh_partial.close()
     fh_sub.close()
     fh_ref.close()
     fh_indel.close()
     fh_mixed.close()
     fh_cutsub.close()
     hist_inframe=dict(hist_inframe)
     hist_frameshift=dict(hist_frameshift)

     return      df_needle_alignment_chunk, effect_vector_insertion,effect_vector_deletion,\
     effect_vector_mutation,effect_vector_any,effect_vector_insertion_mixed,effect_vector_deletion_mixed,\
     effect_vector_mutation_mixed,effect_vector_insertion_hdr,effect_vector_deletion_hdr,effect_vector_mutation_hdr,\
     effect_vector_insertion_noncoding,effect_vector_deletion_noncoding,effect_vector_mutation_noncoding,hist_inframe,\
     hist_frameshift,avg_vector_del_all,avg_vector_ins_all,MODIFIED_FRAMESHIFT,MODIFIED_NON_FRAMESHIFT,NON_INDELS_NON_FRAMESHIFT,\
     SPLICING_SITES_MODIFIED




def add_hist(hist_to_add,hist_global):
    for key,value in hist_to_add.iteritems():
        hist_global[key]+=value
    return hist_global


def slugify(value): #adapted from the Django project

    print("causing error " + value);
    #value = unicodedata.normalize('NFKD', unicode(value)).encode('ascii', 'ignore')
    #value = unicode(re.sub('[-\s]+', '-', value))
    value = unicodedata.normalize('NFKD', value)
    value = str(re.sub('[^\w\s-]', '_', value).strip())
    value = str(re.sub('[^\w\s-]', '_', value).strip())

    return str(value)

def split_paired_end_reads_single_file(fastq_filename,output_filename_r1,output_filename_r2):

    if fastq_filename.endswith('.gz'):
            fastq_handle=gzip.open(fastq_filename)
    else:
            fastq_handle=open(fastq_filename)

    #we cannot use with on gzip with python 2.6 :(
    try:
        fastq_splitted_outfile_r1=gzip.open(output_filename_r1,'w+')
        fastq_splitted_outfile_r2=gzip.open(output_filename_r2,'w+')
        [fastq_splitted_outfile_r1.write(line) if (i % 8 < 4) else fastq_splitted_outfile_r2.write(line) for i, line in enumerate(fastq_handle)]
    except:
        raise Exception('Error handling the splitting operation')

    return output_filename_r1,output_filename_r2


def get_row_around_cut(row,cut_point,offset):
    cut_idx=row['ref_positions'].index(cut_point)
    return  row['Aligned_Sequence'][cut_idx-offset+1:cut_idx+offset+1],row['Reference_Sequence'][cut_idx-offset+1:cut_idx+offset+1],row['UNMODIFIED'],row['%Reads'], row['#Reads']


def get_dataframe_around_cut(df_alleles, cut_point,offset):
    df_alleles_around_cut=pd.DataFrame(list(df_alleles.apply(lambda row: get_row_around_cut(row,cut_point,offset),axis=1).values), columns=['Aligned_Sequence','Reference_Sequence','Unedited','%Reads','#Reads'])
    df_alleles_around_cut=df_alleles_around_cut.groupby(['Aligned_Sequence','Reference_Sequence']).sum().reset_index().set_index('Aligned_Sequence')

    df_alleles_around_cut.sort_values(by='%Reads',inplace=True,ascending=False)
    df_alleles_around_cut['Unedited']=df_alleles_around_cut['Unedited']>0
    return df_alleles_around_cut

#We need to customize the seaborn heatmap class and function
class Custom_HeatMapper(sns.matrix._HeatMapper):

    def __init__(self, data, vmin, vmax, cmap, center, robust, annot, fmt,
                 annot_kws,per_element_annot_kws,cbar, cbar_kws,
                 xticklabels=True, yticklabels=True, mask=None):

        super(Custom_HeatMapper, self).__init__(data, vmin, vmax, cmap, center, robust, annot, fmt,
                 annot_kws, cbar, cbar_kws,
                 xticklabels, yticklabels, mask)


        if annot is not None:
            if per_element_annot_kws is None:
                self.per_element_annot_kws=np.empty_like(annot,dtype=np.object)
                self.per_element_annot_kws[:]=dict()
            else:
                self.per_element_annot_kws=per_element_annot_kws

    #add per element dict to syle the annotatiin
    def _annotate_heatmap(self, ax, mesh):
        """Add textual labels with the value in each cell."""
        mesh.update_scalarmappable()
        xpos, ypos = np.meshgrid(ax.get_xticks(), ax.get_yticks())


        for x, y, m, color, val,per_element_dict  in zip(xpos.flat, ypos.flat,
                                       mesh.get_array(), mesh.get_facecolors(),
                                       self.annot_data.flat,self.per_element_annot_kws.flat):
            #print per_element_dict
            if m is not np.ma.masked:
                l = sns.utils.relative_luminance(color)
                text_color = ".15" if l > .408 else "w"
                annotation = ("{:" + self.fmt + "}").format(val)
                text_kwargs = dict(color=text_color, ha="center", va="center")
                text_kwargs.update(self.annot_kws)
                text_kwargs.update(per_element_dict)

                ax.text(x, y, annotation, **text_kwargs)


    #removed the colobar
    def plot(self, ax, cax, kws):
        """Draw the heatmap on the provided Axes."""
        # Remove all the Axes spines
        sns.utils.despine(ax=ax, left=True, bottom=True)

        # Draw the heatmap
        mesh = ax.pcolormesh(self.plot_data, vmin=self.vmin, vmax=self.vmax,
                             cmap=self.cmap, **kws)

        # Set the axis limits
        ax.set(xlim=(0, self.data.shape[1]), ylim=(0, self.data.shape[0]))

        # Add row and column labels
        ax.set(xticks=self.xticks, yticks=self.yticks)
        xtl = ax.set_xticklabels(self.xticklabels)
        ytl = ax.set_yticklabels(self.yticklabels, rotation="vertical")

        # Possibly rotate them if they overlap
        plt.draw()
        if sns.utils.axis_ticklabels_overlap(xtl):
            plt.setp(xtl, rotation="vertical")
        if sns.utils.axis_ticklabels_overlap(ytl):
            plt.setp(ytl, rotation="horizontal")

        # Add the axis labels
        ax.set(xlabel=self.xlabel, ylabel=self.ylabel)

        # Annotate the cells with the formatted values
        if self.annot:
            self._annotate_heatmap(ax, mesh)




def custom_heatmap(data, vmin=None, vmax=None, cmap=None, center=None, robust=False,
            annot=None, fmt=".2g", annot_kws=None,per_element_annot_kws=None,
            linewidths=0, linecolor="white",
            cbar=True, cbar_kws=None, cbar_ax=None,
            square=False, ax=None, xticklabels=True, yticklabels=True,
            mask=None,
            **kwargs):

    # Initialize the plotter object
    plotter = Custom_HeatMapper(data, vmin, vmax, cmap, center, robust, annot, fmt,
                          annot_kws, per_element_annot_kws,cbar, cbar_kws, xticklabels,
                          yticklabels, mask)

    # Add the pcolormesh kwargs here
    kwargs["linewidths"] = linewidths
    kwargs["edgecolor"] = linecolor

    # Draw the plot and return the Axes
    if ax is None:
        ax = plt.gca()
    if square:
        ax.set_aspect("equal")
    plotter.plot(ax, cbar_ax, kwargs)
    return ax

#def plot_alleles_table(offset_to_plot,reference_seq,cut_point,df_alleles,sgRNA_name,OUTPUT_DIRECTORY,MIN_FREQUENCY=0.0001,MAX_N_ROWS=500):
def plot_alleles_table(offset_to_plot,reference_seq,cut_point,df_alleles,sgRNA_name,OUTPUT_DIRECTORY,MIN_FREQUENCY=0.001,MAX_N_ROWS=150):
    #bp we are plotting on each side SKW
# START FIG 9
    offset_around_cut_to_plot=offset_to_plot

    info('Printing offset of %d' % offset_around_cut_to_plot)

    # make a color map of fixed colors
    alpha=0.5
    get_color=lambda x,y,z: (x/255.0,y/255.0,z/255.0,alpha)
    A_color=get_color(127,201,127)
    T_color=get_color(190,174,212)
    C_color=get_color(253,192,134)
    G_color=get_color(255,255,153)
    INDEL_color=get_color(230,230,230)

    cmap = colors_mpl.ListedColormap([INDEL_color, A_color,T_color,C_color,G_color])

    dna_to_numbers={'-':0,'A':1,'T':2,'C':3,'G':4,'N':5}
    seq_to_numbers= lambda seq: [dna_to_numbers[x] for x in seq]

    X=[]
    annot=[]
    y_labels=[]
    lines=defaultdict(list)
    re_find_indels=re.compile("(-*-)")

    per_element_annot_kws=[]
    idx_row=0
    for idx,row in df_alleles.ix[df_alleles['%Reads']>=MIN_FREQUENCY][:MAX_N_ROWS].iterrows():
        X.append(seq_to_numbers(str.upper(idx)))
        annot.append(list(idx))
        y_labels.append('%.2f%% (%d reads)' % (row['%Reads'],row['#Reads']))

        for p in re_find_indels.finditer(row['Reference_Sequence']):
            lines[idx_row].append((p.start(),p.end()))

        idx_row+=1
        idxs_sub= [i_sub for i_sub in range(len(idx)) if \
                   (row['Reference_Sequence'][i_sub]!=idx[i_sub]) and \
                   (row['Reference_Sequence'][i_sub]!='-') and\
                   (idx[i_sub]!='-')]
        to_append=np.array([{}]*len(idx),dtype=np.object)
        to_append[ idxs_sub]={'weight':'bold', 'color':'black','size':16}
        per_element_annot_kws.append(to_append)

    ref_seq_around_cut=reference_seq[cut_point-offset_around_cut_to_plot+1:cut_point+offset_around_cut_to_plot+1]

    per_element_annot_kws=np.vstack(per_element_annot_kws[::-1])
    ref_seq_hm=np.expand_dims(seq_to_numbers(ref_seq_around_cut),1).T
    ref_seq_annot_hm=np.expand_dims(list(ref_seq_around_cut),1).T
    NEW_SEABORN=np.sum(np.array(map(int,sns.__version__.split('.')))*(100,10,1))>= 80

    if NEW_SEABORN:
        annot=annot[::-1]
        X=X[::-1]

    sns.set_context('poster')

    N_ROWS=len(X)
    N_COLUMNS=offset_around_cut_to_plot*2

    fig=plt.figure(figsize=(offset_around_cut_to_plot*0.6,(N_ROWS+1)*0.6))
    gs1 = gridspec.GridSpec(N_ROWS+1,N_COLUMNS)
    gs2 = gridspec.GridSpec(N_ROWS+1,N_COLUMNS)

    ax_hm_ref=plt.subplot(gs1[0, :])
    ax_hm=plt.subplot(gs2[1:, :])

    custom_heatmap(ref_seq_hm,annot=ref_seq_annot_hm,annot_kws={'size':16},cmap=cmap,fmt='s',ax=ax_hm_ref,vmin=0,vmax=5,square=True)
    custom_heatmap(X,annot=np.array(annot),annot_kws={'size':16},cmap=cmap,fmt='s',ax=ax_hm,square=True, per_element_annot_kws=per_element_annot_kws)

    ax_hm.yaxis.tick_right()
    ax_hm.yaxis.set_ticklabels(y_labels[::-1],rotation=True),
    ax_hm.xaxis.set_ticks([])

    #print lines

    #cut point vertical line
    ax_hm.vlines([offset_around_cut_to_plot],*ax_hm.get_ylim(),linestyles='dashed')

    #create boxes for ins
    for idx,lss in lines.iteritems():
            for ls in lss:
                for l in ls:
                    ax_hm.vlines([l],N_ROWS-idx-1,N_ROWS-idx,color='red',lw=3)

                ax_hm.hlines(N_ROWS-idx-1,ls[0],ls[1],color='red',lw=3)
                ax_hm.hlines(N_ROWS-idx,ls[0],ls[1],color='red',lw=3)

    ax_hm_ref.yaxis.tick_right()
    ax_hm_ref.xaxis.set_ticks([])
    ax_hm_ref.yaxis.set_ticklabels(['Reference'],rotation=True)

    gs2.update(left=0,right=1, hspace=0.05,wspace=0,top=1*(((N_ROWS)*1.13))/(N_ROWS))
    gs1.update(left=0,right=1, hspace=0.05,wspace=0,)

    sns.set_context(rc={'lines.markeredgewidth': 1,'mathtext.fontset' : 'stix','text.usetex':True,'text.latex.unicode':True} )

    proxies = [matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='black',
                    mec='none', marker=r'$\mathbf{{{}}}$'.format('bold'),ms=18),
               matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                    mec='red', marker='s',ms=8,markeredgewidth=2.5),
              matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                    mec='black', marker='_',ms=2,),
              matplotlib.lines.Line2D([0], [1], linestyle='--',c='black',ms=6)] #
    descriptions=['Substitutions','Insertions','Deletions','Predicted cleavage position']
    ax_hm_ref.legend(proxies, descriptions, numpoints=1, markerscale=2, loc='center', bbox_to_anchor=(0.5, 4),ncol=1)

    _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename)

    plt.savefig(_jp('9.Alleles_around_cut_site_for_%s.pdf' % args.name),bbox_inches='tight')
    if args.save_also_png:
        plt.savefig(_jp('9.Alleles_around_cut_site_for_%s.png' % args.name),bbox_inches='tight',pad=1)
    print('Plot 9 Done\n')
# END FIG 9


def main():
    try:
             print('  \n~~~cortado~~~')
             print('-Analysis of CRISPR/Cas9 outcomes from deep sequencing data-')
             print('Reimplemented from CRISPResso Version 1.0.8\n')
             print('''
                      ))
                     ((
                      ))
                    _____
                 C\|~~~~~|
                    \___/
             ''')
             #print'\n[Luca Pinello 2015, send bugs, suggestions or *green coffee* to lucapinello AT gmail DOT com]\n\n',

             #global variables for the multiprocessing

             global args
             global include_idxs
             global len_amplicon
             global exon_positions
             global global_sub_count
             global partial_correction
             global seq_error
             global splicing_positions
             global hdr_vector 
             global fh_ref,fh_sub,fh_hdr,fh_index,fh_mixed,fh_cutsub,fh_indel,fh_partial
             global ref_donor_diffs
            
             global_sub_count = 0
             seq_error = 0 
             partial_correction = 0 
             parser = argparse.ArgumentParser(description='cortado parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
             parser.add_argument('-r1','--fastq_r1', type=str,  help='First fastq file', required=True,default='Fastq filename' )
             parser.add_argument('-r2','--fastq_r2', type=str,  help='Second fastq file for paired end reads',default='')
             parser.add_argument('-a','--amplicon_seq', type=str,  help='Amplicon Sequence', required=True)

             #optional
             parser.add_argument('-g','--guide_seq',  help="sgRNA sequence, if more than one, please separate by comma/s. Note that the sgRNA needs to be input as the guide RNA sequence (usually 20 nt) immediately adjacent to but not including the PAM sequence (5' of NGG for SpCas9). If the PAM is found on the opposite strand with respect to the Amplicon Sequence, ensure the sgRNA sequence is also found on the opposite strand. The cortado convention is to depict the expected cleavage position using the value of the parameter cleavage_offset nt  3' from the end of the guide. In addition, the use of alternate nucleases to SpCas9 is supported. For example, if using the Cpf1 system, enter the sequence (usually 20 nt) immediately 3' of the PAM sequence and explicitly set the cleavage_offset parameter to 1, since the default setting of -3 is suitable only for SpCas9.", default='')
             parser.add_argument('-e','--expected_hdr_amplicon_seq',  help='Amplicon sequence expected after hdr', default='')
             parser.add_argument('-d','--donor_seq',  help='Donor Sequence. This optional input comprises a subsequence of the expected hdr amplicon to be highlighted in plots.', default='')
             parser.add_argument('-c','--coding_seq',  help='Subsequence/s of the amplicon sequence covering one or more coding sequences for the frameshift analysis.If more than one (for example, split by intron/s), please separate by comma.', default='')
             parser.add_argument('-q','--min_average_read_quality', type=int, help='Minimum average quality score (phred33) to keep a read', default=0)
             parser.add_argument('-s','--min_single_bp_quality', type=int, help='Minimum single bp score (phred33) to keep a read', default=0)
             parser.add_argument('--min_identity_score', type=float, help='Minimum identity score for the alignment', default=60.0)
             parser.add_argument('-n','--name',  help='Output name', default='')
             parser.add_argument('-k','--skip_aln',  help='Skip alignment step', default='')
             parser.add_argument('--main_site',  type=int, help='Main HDR site (numbered left to right)', default=1)
             parser.add_argument('--all_edits',  help='Output column with number of reads with all edits', action='store_true')
             parser.add_argument('-o','--output_folder',  help='', default='')
             parser.add_argument('--split_paired_end',help='Splits a single fastq file contating paired end reads in two files before running CRISPResso',action='store_true')
             parser.add_argument('--trim_sequences',help='Enable the trimming of Illumina adapters with Trimmomatic',action='store_true')
             parser.add_argument('--trimmomatic_options_string', type=str, help='Override options for Trimmomatic',default=' ILLUMINACLIP:%s:2:30:10' % get_data('adapters.fa'))
             #parser.add_argument('--trimmomatic_options_string', type=str, help='Override options for Trimmomatic',default=' ILLUMINACLIP:%s:0:90:10:0:true MINLEN:40' % get_data('NexteraPE-PE.fa'))
             parser.add_argument('--min_paired_end_reads_overlap',  type=int, help='Minimum required overlap length between two reads to provide a confident overlap. ', default=75)
             parser.add_argument('--hide_mutations_outside_window_NHEJ',help='This parameter allows to visualize only the mutations overlapping the cleavage site and used to classify a read as NHEJ. This parameter has no effect on the quanitification of the NHEJ. It  may be helpful to mask a pre-existing and known mutations or sequencing errors outside the window used for quantification of NHEJ events.',action='store_true')
             parser.add_argument('-w','--window_around_sgrna', type=int, help='Window(s) in bp around the cleavage position (half on on each side) as determined by the provide guide RNA sequence to quantify the indels. Any indels outside this window are excluded. A value of 0 disables this filter.', default=1)
             parser.add_argument('--cleavage_offset', type=int, help="Cleavage offset to use within respect to the 3' end of the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. The default is -3 and is suitable for the SpCas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 this parameter would be set to 1.", default=-3)
             parser.add_argument('--exclude_bp_from_left', type=int, help='Exclude bp from the left side of the amplicon sequence for the quantification of the indels', default=15)
             parser.add_argument('--exclude_bp_from_right', type=int, help='Exclude bp from the right side of the amplicon sequence for the quantification of the indels', default=15)
             parser.add_argument('--hdr_perfect_alignment_threshold',  type=float, help='Sequence homology %% for an HDR occurrence', default=98.0)
             parser.add_argument('--ignore_substitutions',help='Ignore substitutions events for the quantification and visualization',action='store_true')
             parser.add_argument('--ignore_insertions',help='Ignore insertions events for the quantification and visualization',action='store_true')
             parser.add_argument('--ignore_deletions',help='Ignore deletions events for the quantification and visualization',action='store_true')
             parser.add_argument('--needle_options_string',type=str,help='Override options for the Needle aligner',default='-gapopen=10 -gapextend=0.5  -awidth3=5000')
             parser.add_argument('--keep_intermediate',help='Keep all the  intermediate files',action='store_true')
             parser.add_argument('--dump',help='Dump numpy arrays and pandas dataframes to file for debugging purposes',action='store_true')
             parser.add_argument('--save_also_png',help='Save also .png images additionally to .pdf files',action='store_true')
             parser.add_argument('-p','--n_processes',type=int, help='Specify the number of processes to use for the quantification.\
             Please use with caution since increasing this parameter will increase significantly the memory required to run CRISPResso.',default=1)
             parser.add_argument('--offset_around_cut_to_plot',  type=int, help='Offset to use to summarize alleles around the cut site in the alleles table plot.', default=20)
             parser.add_argument('--min_frequency_alleles_around_cut_to_plot', type=float, help='Minimum %% reads required to report an allele in the alleles table plot.', default=0.01)
             parser.add_argument('--max_rows_alleles_around_cut_to_plot',  type=int, help='Maximum number of rows to report in the alleles table plot. ', default=250)
             args = parser.parse_args()

             # vector of subs is zero based.
             args.main_site = args.main_site - 1
             #check files
             check_file(args.fastq_r1)
             if args.fastq_r2:
                     check_file(args.fastq_r2)

             #normalize name and remove not allowed characters
             if args.name:
                 clean_name=slugify(args.name)
                 if args.name!= clean_name:
                        warn('The specified name %s contained characters not allowed and was changed to: %s' % (args.name,clean_name))
                        args.name=clean_name

             #amplicon sequence check
             #make evetything uppercase!
             args.amplicon_seq=args.amplicon_seq.strip().upper()
             wrong_nt=find_wrong_nt(args.amplicon_seq)
             if wrong_nt:
                 raise NTException('The amplicon sequence contains wrong characters:%s' % ' '.join(wrong_nt))

             len_amplicon=len(args.amplicon_seq)
             #hdr_vector=np.zeros(5)
             hdr_vector=np.zeros(50) #HERE is what I changed for mohan's data

             if args.guide_seq:
                     cut_points=[]
                     sgRNA_intervals=[]
                     offset_plots=[]
                     sgRNA_sequences=[]

                     args.guide_seq=args.guide_seq.strip().upper()
                     for current_guide_seq in args.guide_seq.split(','):
                         if current_guide_seq in args.amplicon_seq:
                            offset_plots.append(1)
                         else:
                            offset_plots.append(0)

                         wrong_nt=find_wrong_nt(current_guide_seq)
                         if wrong_nt:
                            raise NTException('The sgRNA sequence contains wrong characters:%s'  % ' '.join(wrong_nt))

                         offset_fw=args.cleavage_offset+len(current_guide_seq)-1
                         offset_rc=(-args.cleavage_offset)-1
                         cut_points+=[m.start() + offset_fw for m in re.finditer(current_guide_seq, args.amplicon_seq)]+[m.start() + offset_rc for m in re.finditer(reverse_complement(current_guide_seq), args.amplicon_seq)]
                         sgRNA_intervals+=[(m.start(),m.start()+len(current_guide_seq)-1) for m in re.finditer(current_guide_seq, args.amplicon_seq)]+[(m.start(),m.start()+len(current_guide_seq)-1) for m in re.finditer(reverse_complement(current_guide_seq), args.amplicon_seq)]
                         sgRNA_sequences.append(current_guide_seq)

                     offset_plots=np.array(offset_plots)

                     if not cut_points:
                         raise SgRNASequenceException('The guide sequence/s provided is(are) not present in the amplicon sequence! \n\nPlease check your input!')
                     else:
                         info('Cut site from guide seq:%s' % cut_points)

             else:
                     cut_points=[]
                     sgRNA_intervals=[]
                     offset_plots=np.array([])
                     sgRNA_sequences=[]

             # if cut_points > 1/2 len of amplicon, then offset = amplen - cutsite, else cutsite
             #args.offset_around_cut_to_plot = cut_points[0] - 1


             if args.expected_hdr_amplicon_seq:
                     args.expected_hdr_amplicon_seq=args.expected_hdr_amplicon_seq.strip().upper()
                     #ref_donor_diffs=[i for i in xrange(len(args.expected_hdr_amplicon_seq)) if args.expected_hdr_amplicon_seq[i] != args.amplicon_seq[i]]

                     if args.expected_hdr_amplicon_seq == args.amplicon_seq:
                         raise AmpliconEqualDonorException('The amplicon sequence expected after an hdr and the reference amplicon cannot be the same! \n\nPlease check your input!')

                     wrong_nt=find_wrong_nt(args.expected_hdr_amplicon_seq)
                     if wrong_nt:
                        raise NTException('The amplicon sequence expected after an hdr contains wrong characters:%s' % ' '.join(wrong_nt))

                     #if len(args.expected_hdr_amplicon_seq)!=len(args.amplicon_seq):
                     aligned_ref,aligned_exp=pairwise2.align.globalxx (args.amplicon_seq,args.expected_hdr_amplicon_seq)[0][:2]
                     identity_ref_rep=sum([1.0 for a,b in zip(aligned_ref,aligned_exp)  if a==b  ])/len(aligned_ref)*100


                     if identity_ref_rep < args.min_identity_score:
                         raise DonorSequenceException('The amplicon sequence expected after an hdr should be provided as the reference amplicon sequence with the relevant part of the donor sequence replaced, and not just as the donor sequence. \n\nPlease check your input!')

             if args.donor_seq:
                     args.donor_seq=args.donor_seq.strip().upper()
                     wrong_nt=find_wrong_nt(args.donor_seq)
                     if wrong_nt:
                         raise NTException('The donor sequence contains wrong characters:%s' % ' '.join(wrong_nt))

                     if args.donor_seq not in args.expected_hdr_amplicon_seq:
                         raise CoreDonorSequenceNotContainedException('The donor sequence provided is not present in the expected hdr amplicon sequence, or the expected hdr amplicon sequence parameter (-e) is not defined.  \n\nPlease check your input!')

                     positions_core_donor_seq=[(m.start(),m.start()+len(args.donor_seq)) for m in re.finditer('(?=%s)' % args.donor_seq, args.expected_hdr_amplicon_seq)]
                     if len(positions_core_donor_seq)>1:
                         raise CoreDonorSequenceNotUniqueException('The donor sequence provided is not unique in the expected hdr amplicon sequence.  \n\nPlease check your input!')
                     core_donor_seq_st_en=positions_core_donor_seq[0]



             ###FRAMESHIFT SUPPORT###
             if args.coding_seq:
                    PERFORM_FRAMESHIFT_ANALYSIS=True
                    exon_positions=set()
                    exon_intervals=[]
                    splicing_positions=[]

                    for exon_seq in args.coding_seq.strip().upper().split(','):
                        #check for wrong NT
                        wrong_nt=find_wrong_nt(exon_seq)
                        if wrong_nt:
                            raise NTException('The coding sequence contains wrong characters:%s' % ' '.join(wrong_nt))

                        st_exon=args.amplicon_seq.find(exon_seq)
                        if  st_exon<0:
                            raise ExonSequenceException('The coding subsequence/s provided:%s is(are) not contained in the amplicon sequence.' % exon_seq)
                        en_exon=st_exon+len(exon_seq) #this do not include the upper bound as usual in python
                        exon_intervals.append((st_exon,en_exon))
                        exon_positions=exon_positions.union(set(range(st_exon,en_exon)))
                        #consider 2 base pairs before and after each exon
                        splicing_positions+=[max(0,st_exon-2),max(0,st_exon-1),min(len_amplicon-1, en_exon),min(len_amplicon-1, en_exon+1)]

                    exon_positions=sorted(exon_positions)
                    #protect from the wrong splitting of exons by the users to avoid false splicing sites
                    splicing_positions=set(splicing_positions).difference(exon_positions)
             else:
                    PERFORM_FRAMESHIFT_ANALYSIS=False


             #we have insertions/deletions that change the concatenated exon sequence lenght and the difference between the final sequence
             #and the original sequence lenght is not a multiple of 3
             MODIFIED_FRAMESHIFT=0
             #we have insertions/deletions that change the concatenated exon sequence lenght and the difference between the final sequence
             #and the original sequence lenght is a multiple of 3. We are in this case also when no indels are present but we have
             #substitutions
             MODIFIED_NON_FRAMESHIFT=0
             #we don't touch the exons at all, the read can be still modified tough..
             NON_INDELS_NON_FRAMESHIFT=0
             SPLICING_SITES_MODIFIED=0

             ################

             get_name_from_fasta=lambda  x: os.path.basename(x).replace('.fastq','').replace('.gz','')

             if not args.name:
                     if args.fastq_r2!='':
                             database_id='%s_%s' % (get_name_from_fasta(args.fastq_r1),get_name_from_fasta(args.fastq_r2))
                     else:
                             database_id='%s' % get_name_from_fasta(args.fastq_r1)
             else:
                     database_id=args.name

             OUTPUT_DIRECTORY='%s' % database_id

             if args.output_folder:
                     OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)

             _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory
             log_filename=_jp('cortado_running_log.txt')

             try:
                     os.makedirs(OUTPUT_DIRECTORY)
                     info('Creating Folder %s' % OUTPUT_DIRECTORY)
             except:
                     warn('Folder %s already exists, overwriting old output' % OUTPUT_DIRECTORY)

             finally:
                     logging.getLogger().addHandler(logging.FileHandler(log_filename))
                     with open(log_filename,'w+') as outfile:
                         outfile.write('Cut site from guide seq:%s\n' % cut_points)
                         outfile.write('[Command used]:\ncortado %s\n\n[Execution log]:\n' % ' '.join(sys.argv))

             filename=os.path.join(os.path.abspath(args.output_folder),args.name,'hdr_reads.txt')
             fh_hdr = open(filename,"w")
             filename=os.path.join(os.path.abspath(args.output_folder),args.name,'partial_reads.txt')
             fh_partial = open(filename,"w")
             filename=os.path.join(os.path.abspath(args.output_folder),args.name,'indel_reads.txt')
             fh_indel = open(filename,"w")
             filename=os.path.join(os.path.abspath(args.output_folder),args.name,'mixed_reads.txt')
             fh_mixed = open(filename,"w")
             filename=os.path.join(os.path.abspath(args.output_folder),args.name,'ref_reads.txt')
             fh_ref = open(filename,"w")
             filename=os.path.join(os.path.abspath(args.output_folder),args.name,'cutsub_reads.txt')
             fh_cutsub = open(filename,"w")
             filename=os.path.join(os.path.abspath(args.output_folder),args.name,'sub_reads.txt')
             fh_sub = open(filename,"w")

             if args.split_paired_end:
                if args.fastq_r2!='':
                    raise Exception('The option --split_paired_end is available only when a single fastq file is specified!')
                else:
                    info('Splitting paired end single fastq file in two files...')
                    args.fastq_r1,args.fastq_r2=split_paired_end_reads_single_file(args.fastq_r1,
                        output_filename_r1=_jp(os.path.basename(args.fastq_r1.replace('.fastq','')).replace('.gz','')+'_splitted_r1.fastq.gz'),
                        output_filename_r2=_jp(os.path.basename(args.fastq_r1.replace('.fastq','')).replace('.gz','')+'_splitted_r2.fastq.gz'),)
                    splitted_files_to_remove=[args.fastq_r1,args.fastq_r2]
                    print('==>Done.')

             if args.min_average_read_quality>0 or args.min_single_bp_quality>0:
                info('Filtering reads with average bp quality < %d and single bp quality < %d ...' % (args.min_average_read_quality,args.min_single_bp_quality))
                if args.fastq_r2!='':
                        args.fastq_r1,args.fastq_r2=filter_pe_fastq_by_qual(args.fastq_r1,
                             args.fastq_r2,
                             output_filename_r1=_jp(os.path.basename(args.fastq_r1.replace('.fastq','')).replace('.gz','')+'_filtered.fastq.gz'),
                             output_filename_r2=_jp(os.path.basename(args.fastq_r2.replace('.fastq','')).replace('.gz','')+'_filtered.fastq.gz'),
                             min_bp_quality=args.min_average_read_quality,
                             min_single_bp_quality=args.min_single_bp_quality,
                                                                         )
                else:
                        args.fastq_r1=filter_se_fastq_by_qual(args.fastq_r1,
                             output_filename=_jp(os.path.basename(args.fastq_r1).replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz'),
                             min_bp_quality=args.min_average_read_quality,
                             min_single_bp_quality=args.min_single_bp_quality,
                                                                   )

             if not args.skip_aln:
                 if args.fastq_r2=='': #single end reads
                     #check if we need to trim
                     if not args.trim_sequences:
                         #create a symbolic link
                         symlink_filename=_jp(os.path.basename(args.fastq_r1))
                         force_symlink(os.path.abspath(args.fastq_r1),symlink_filename)
                         output_forward_filename=symlink_filename
                     else:
                         output_forward_filename=_jp('reads.trimmed.fq.gz')
                         #Trimming with trimmomatic
                         cmd='java -jar %s SE -threads 2 -phred33 %s  %s %s >>%s 2>&1'\
                         % (get_data('trimmomatic-0.33.jar'),args.fastq_r1,
                            output_forward_filename,
                            args.trimmomatic_options_string.replace('NexteraPE-PE.fa','TruSeq3-SE.fa'),
                            log_filename)
                         TRIMMOMATIC_STATUS=sb.call(cmd,shell=True)
    
                         if TRIMMOMATIC_STATUS:
                                 raise TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')
    
                     processed_output_filename=output_forward_filename
    
                 else:#paired end reads case
    
                     if not args.trim_sequences:
                         output_forward_paired_filename=args.fastq_r1
                         output_reverse_paired_filename=args.fastq_r2
                     else:
                         print('Trimming sequences with Trimmomatic...')
                         output_forward_paired_filename=_jp('output_forward_paired.fq.gz')
                         output_forward_unpaired_filename=_jp('output_forward_unpaired.fq.gz')
                         output_reverse_paired_filename=_jp('output_reverse_paired.fq.gz')
                         output_reverse_unpaired_filename=_jp('output_reverse_unpaired.fq.gz')
    
                         #Trimming with trimmomatic
                         cmd='java -jar %s PE -threads 2 -phred33 %s  %s %s  %s  %s  %s %s >>%s 2>&1'\
                         % (get_data('trimmomatic-0.33.jar'),
                                 args.fastq_r1,args.fastq_r2,output_forward_paired_filename,
                                 output_forward_unpaired_filename,output_reverse_paired_filename,
                                 output_reverse_unpaired_filename,args.trimmomatic_options_string,log_filename)
                         info(cmd)
                         TRIMMOMATIC_STATUS=sb.call(cmd,shell=True)
                         if TRIMMOMATIC_STATUS:
                                 raise TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')
    
                         print('==>Done.')
    
                     print ('Estimating average read length...')
                     if get_n_reads_fastq(output_forward_paired_filename):
                         avg_read_length=get_avg_read_lenght_fastq(output_forward_paired_filename)
                         std_fragment_length=int(len_amplicon*0.1)
                     else:
                        raise NoReadsAfterQualityFiltering('No reads survived the average or single bp quality filtering.')
    
                     #Merging with Flash
                     print('Merging paired sequences with Flash...')
                     cmd='flash %s %s -t 2 --allow-outies --max-overlap 300 --min-overlap 10 -f %d -r %d -s %d  -z -d %s >>%s 2>&1' %\
                     (output_forward_paired_filename,
                      output_reverse_paired_filename,
                      len_amplicon,avg_read_length,
                      std_fragment_length,
                      OUTPUT_DIRECTORY,log_filename)

                     info(cmd)
                     FLASH_STATUS=sb.call(cmd,shell=True)
                     if FLASH_STATUS:
                         raise FlashException('Flash failed to run, please check the log file.')
    
                     print('==>Done joining reads.')
                     flash_hist_filename=_jp('out.hist')
                     flash_histogram_filename=_jp('out.histogram')
                     flash_not_combined_1_filename=_jp('out.notCombined_1.fastq.gz')
                     flash_not_combined_2_filename=_jp('out.notCombined_2.fastq.gz')
    
                     processed_output_filename=_jp('out.extendedFrags.fastq.gz')
                     #================= reorient =========================================
                     # Reorient all the reads in the same dir
                     oriented_output_filename=_jp('out.oriented.fastq.gz')
                     cmd='reorient.pl %s %s %s' %\
                     (processed_output_filename,
                      args.amplicon_seq,
                      oriented_output_filename)
                     info(cmd)
                     REORIENT_STATUS=sb.call(cmd,shell=True)
                     if REORIENT_STATUS:
                         raise ReorientException('Reorient failed to run, please check the log file.')
                     processed_output_filename=oriented_output_filename
                     #exit()
                     #================= reorient =========================================
                 #count reads
                 N_READS_INPUT=get_n_reads_fastq(args.fastq_r1)
                 N_READS_AFTER_PREPROCESSING=get_n_reads_fastq(processed_output_filename)
                 print("processed output filesnmae: %s %d" % (processed_output_filename,N_READS_AFTER_PREPROCESSING))
                 if N_READS_AFTER_PREPROCESSING < 1000:			#SKW less than 1K, no sensitivity
                     err_str = 'NOT ENOUGH READS: {0} in fastq, {1} after reorienting'.format(N_READS_INPUT,N_READS_AFTER_PREPROCESSING)
                     raise NoReadsAfterQualityFiltering(err_str)
                 #write statistics
                 #with open(_jp('Mapping_statistics.txt'),'w+') as outfile:
                 #    outfile.write('READS IN INPUTS:%d\nREADS AFTER PREPROCESSING:%d\nREADS ALIGNED:%d' % (N_READS_INPUT,N_READS_AFTER_PREPROCESSING,N_ALIGNED))

             print('Preparing files for the alignment...')
             #parsing flash output and prepare the files for alignment

             needle_output_filename=_jp('needle_output_%s.txt.gz' % database_id)
             needle_output_repair_filename=_jp('needle_output_repair_%s.txt.gz' % database_id)

             #write .fa file only for amplicon the rest we pipe trough awk on the fly!
             if not args.skip_aln:
                 database_fasta_filename=_jp('%s_database.fa' % database_id)
                 with open(database_fasta_filename,'w+') as outfile:
                         outfile.write('>%s\n%s\n' % (database_id,args.amplicon_seq))

                 if args.expected_hdr_amplicon_seq:
                     database_repair_fasta_filename=_jp('%s_database_repair.fa' % database_id)
                     with open(database_repair_fasta_filename,'w+') as outfile:
                             outfile.write('>%s\n%s\n' % (database_id,args.expected_hdr_amplicon_seq))
             print('==>Done.')

             def parse_needle_output(needle_filename,name='seq',just_score=False):
                 needle_data=[]

                 try:
                     needle_infile = gzip.open(needle_filename,'rt')
                     line = needle_infile.readline()
                     while line:
                         while line and ('# Aligned_sequences' not  in line):
                             line=needle_infile.readline()

                         if line:
                                     #print line
                                     needle_infile.readline() #skip another line
                                     line=needle_infile.readline()
                                     id_seq=line.split()[-1].replace('_',':')
                                     for _ in range(5):
                                             needle_infile.readline()
                                     line=needle_infile.readline()
                                     identity_seq=eval(line.strip().split(' ')[-1].replace('%','').replace(')','').replace('(',''))
                                     for _ in range(7):
                                         needle_infile.readline()

                                     line=needle_infile.readline()
                                     aln_ref_seq=line.split()[2]

                                     aln_str=needle_infile.readline()[21:].rstrip('\n')
                                     line=needle_infile.readline()
                                     aln_query_seq=line.split()[2]
                                     aln_query_len=line.split()[3]
                                     if just_score:
                                             needle_data.append([id_seq,identity_seq,aln_ref_seq])
                                     else:
                                             needle_data.append([id_seq,identity_seq,aln_query_len,aln_ref_seq,aln_str,aln_query_seq])

                     if just_score:
                         needle_infile.close()
                         return pd.DataFrame(needle_data,columns=['ID','score_'+name,'aln_donor_seq']).set_index('ID')
                     else:
                         needle_infile.close()
                         return pd.DataFrame(needle_data,columns=['ID','score_'+name,'length','ref_seq','align_str','align_seq']).set_index('ID')
                 except:
                     raise NeedleException('Failed to parse the output of needle!')

             print('Aligning sequences with needle...')

             if not args.skip_aln:
                 cmd=(('cat %s |'% processed_output_filename )+\
                 (' gunzip |' if processed_output_filename.endswith('.gz') else ' '))+\
                 r''' awk 'NR % 4 == 1 {print ">" $0} NR % 4 ==2 {print $0}' '''+\
                 " | sed 's/:/_/g' | needle -asequence=%s -bsequence=/dev/stdin -outfile=/dev/stdout %s 2>> %s  | gzip >%s"\
                 %(database_fasta_filename,args.needle_options_string,log_filename,needle_output_filename)

                 NEEDLE_OUTPUT=sb.call(cmd,shell=True)
                 if NEEDLE_OUTPUT:
                     raise NeedleException('Needle failed to run, please check the log file.')

             #If we have a donor sequence we just compare the fq in the two cases and see which one alignes better
             if args.expected_hdr_amplicon_seq:
                 if not args.skip_aln:
                     cmd_repair=(('cat %s |'% processed_output_filename )+\
                     (' gunzip |' if processed_output_filename.endswith('.gz') else ' '))+\
                     r''' awk 'NR % 4 == 1 {print ">" $0} NR % 4 ==2 {print $0}' '''+\
                     " | sed 's/:/_/g' | needle -asequence=%s -bsequence=/dev/stdin -outfile=/dev/stdout %s 2>> %s  | gzip >%s"\
                     %(database_repair_fasta_filename,args.needle_options_string,log_filename,needle_output_repair_filename)

                     NEEDLE_OUTPUT=sb.call(cmd_repair,shell=True)
                     if NEEDLE_OUTPUT:
                             raise NeedleException('Needle failed to run, please check the log file.')
                     print('==>Done.')
 
             # START HERE IF SKIPPING ALIGNMENT
             #merge the flow
             if args.expected_hdr_amplicon_seq:
                    df_database=parse_needle_output(needle_output_filename,'ref')
                    df_database_repair=parse_needle_output(needle_output_repair_filename,'repaired',just_score=True)
                    df_database_and_repair=df_database.join(df_database_repair)

                    del df_database
                    del df_database_repair

                    #filter bad alignments

                    N_TOTAL_ALSO_UNALIGNED=df_database_and_repair.shape[0]*1.0

                    #find reads that failed to align and try on the reverse complement
                    #sr_not_aligned=df_database_and_repair.ix[(df_database_and_repair.score_ref <args.min_identity_score)\
                    sr_not_aligned=df_database_and_repair.loc[(df_database_and_repair.score_ref <args.min_identity_score)\
                                      & (df_database_and_repair.score_ref< args.min_identity_score)]\
                                     .align_seq.apply(lambda x: x.replace('_',''))

                    print("did we get HERE???????????")
                    #filter out not aligned reads
                    df_database_and_repair=\
                    df_database_and_repair.ix[\
                        (df_database_and_repair.score_ref>args.min_identity_score)\
                        |(df_database_and_repair.score_repaired>args.min_identity_score)]

                    df_database_and_repair['score_diff']=df_database_and_repair.score_ref-df_database_and_repair.score_repaired
                    df_needle_alignment=df_database_and_repair
                    del df_database_and_repair

             else:
                    df_needle_alignment=parse_needle_output(needle_output_filename,'ref')
                    N_TOTAL_ALSO_UNALIGNED=df_needle_alignment.shape[0]*1.0

                    sr_not_aligned=df_needle_alignment.ix[(df_needle_alignment.score_ref <args.min_identity_score)]\
                                     .align_seq.apply(lambda x: x.replace('_',''))
                    #filter out not aligned reads
                    df_needle_alignment=df_needle_alignment.ix[df_needle_alignment.score_ref>args.min_identity_score]


             #check if the not aligned reads are in the reverse complement
             if sr_not_aligned.count():
                 #write fastq_not_aligned
                 fasta_not_aligned_filename=_jp('not_aligned_amplicon_forward.fa.gz')

                 if not args.skip_aln:
                     outfile=gzip.open(fasta_not_aligned_filename,'w+')
                     for x in sr_not_aligned.iteritems():
                        outfile.write('>%s\n%s\n' % (x[0],x[1]))

                 #write reverse complement of ampl and expected amplicon
                 database_rc_fasta_filename=_jp('%s_database_rc.fa' % database_id)
                 needle_output_rc_filename=_jp('needle_output_rc_%s.txt.gz' % database_id)

                 print('Align sequences to reverse complement of the amplicon...')

                 if not args.skip_aln:
                     with open(database_rc_fasta_filename,'w+') as outfile:
                         outfile.write('>%s\n%s\n' % (database_id,reverse_complement(args.amplicon_seq)))

                 if args.expected_hdr_amplicon_seq:
                         database_repair_rc_fasta_filename=_jp('%s_database_repair_rc.fa' % database_id)
                         needle_output_repair_rc_filename=_jp('needle_output_repair_rc_%s.txt.gz' % database_id)

                         if not args.skip_aln:
                             with open(database_repair_rc_fasta_filename,'w+') as outfile:
                                 outfile.write('>%s\n%s\n' % (database_id,reverse_complement(args.expected_hdr_amplicon_seq)))
                 print('==>Done.')

                 #Now we do the alignment
                 if not args.skip_aln:
                     cmd="zcat < %s | sed 's/:/_/g' | needle -asequence=%s -bsequence=/dev/stdin -outfile=/dev/stdout %s 2>> %s  | gzip >%s"\
                     %(fasta_not_aligned_filename,database_rc_fasta_filename,args.needle_options_string,log_filename,needle_output_rc_filename)

                     NEEDLE_OUTPUT=sb.call(cmd,shell=True)
                     if NEEDLE_OUTPUT:
                         raise NeedleException('Needle failed to run, please check the log file.')

                 if args.expected_hdr_amplicon_seq:
                    if not args.skip_aln:
                        cmd="zcat < %s | sed 's/:/_/g' | needle -asequence=%s -bsequence=/dev/stdin -outfile=/dev/stdout %s 2>> %s  | gzip >%s"\
                        %(fasta_not_aligned_filename,database_repair_rc_fasta_filename,args.needle_options_string,log_filename,needle_output_repair_rc_filename)

                        NEEDLE_OUTPUT=sb.call(cmd,shell=True)
                        if NEEDLE_OUTPUT:
                           raise NeedleException('Needle failed to run, please check the log file.')

                 #merge the flow rev
                 if args.expected_hdr_amplicon_seq:
                            df_database_rc=parse_needle_output(needle_output_rc_filename,'ref')
                            df_database_repair_rc=parse_needle_output(needle_output_repair_rc_filename,'repaired',just_score=True)
                            df_database_and_repair_rc=df_database_rc.join(df_database_repair_rc)

                            del df_database_rc
                            del df_database_repair_rc

                            #filter bad alignments also to rc
                            df_database_and_repair_rc=\
                            df_database_and_repair_rc.ix[\
                                (df_database_and_repair_rc.score_ref>args.min_identity_score)\
                                |(df_database_and_repair_rc.score_repaired>args.min_identity_score)]

                            df_database_and_repair_rc['score_diff']=df_database_and_repair_rc.score_ref-df_database_and_repair_rc.score_repaired
                            df_needle_alignment_rc=df_database_and_repair_rc
                            del df_database_and_repair_rc

                 else:
                            df_needle_alignment_rc=parse_needle_output(needle_output_rc_filename,'ref')

                            #filter out not aligned reads
                            df_needle_alignment_rc=df_needle_alignment_rc.ix[df_needle_alignment_rc.score_ref>args.min_identity_score]

                 #reverse complement and invert the align string so we have everything in the positive strand
                 df_needle_alignment_rc['ref_seq']=df_needle_alignment_rc['ref_seq'].apply(reverse_complement)
                 df_needle_alignment_rc['align_seq']=df_needle_alignment_rc['align_seq'].apply(reverse_complement)
                 df_needle_alignment_rc['align_str']=df_needle_alignment_rc['align_str'].apply(lambda x: x[::-1])


                 #fix for duplicates when rc alignment
                 df_needle_alignment_rc.index=map(lambda x:'_'.join([x,'RC']),df_needle_alignment_rc.index)

                 #append the RC reads to the aligned reads in the original orientation
                 df_needle_alignment=df_needle_alignment.append(df_needle_alignment_rc)

                 del df_needle_alignment_rc

             #check for duplicates
             try:
                assert df_needle_alignment.shape[0]== df_needle_alignment.index.unique().shape[0]
             except:
                raise DuplicateSequenceIdException('The .fastq file/s contain/s duplicate sequence IDs')


             #Initializations
             print('Quantifying indels/substitutions...')
             df_needle_alignment['UNMODIFIED']=(df_needle_alignment.score_ref==100)

             #the rest we have to look one by one to potentially exclude regions
             df_needle_alignment['MIXED']=False
             df_needle_alignment['HDR']=False
             df_needle_alignment['PARTIAL']=False
             df_needle_alignment['CUTSUB']=False
             df_needle_alignment['NHEJ']=False

             df_needle_alignment['n_mutated']=0
             df_needle_alignment['n_inserted']=0
             df_needle_alignment['n_deleted']=0
             df_needle_alignment['n_substituted']=0

             N_ALIGNED=df_needle_alignment.shape[0]*1.0

             if N_ALIGNED < 10:
                 raise NoReadsAlignedException('Almost no sequences aligned, please check your amplicon sequence')
                 error('No sequences aligned')

             #remove the mutations in bp equal to 'N'
             if 'N' in args.amplicon_seq:
               info('Your amplicon sequence contains one or more N, excluding these bp for the indel quantification...')
               def ignore_N_in_alignment(row):
                   row['align_str']=''.join([('|' if  (row['ref_seq'][idx]=='N') else c ) for idx,c in enumerate(row['align_str'])])
                   if len(set(row['align_str']))==1:
                       row['UNMODIFIED']=True
                   return row
               df_needle_alignment=df_needle_alignment.apply(ignore_N_in_alignment,axis=1)



             #####QUANTIFICATION START
             def compute_ref_positions(ref_seq):
                 pos_idxs=[]
                 idx=0
                 for c in ref_seq:
                     if c in set(['A','T','C','G','N']):
                         pos_idxs.append(idx)
                         idx+=1
                     else:
                         if idx==0:
                             pos_idxs.append(-1)
                         else:
                             pos_idxs.append(-idx)
                 return np.array(pos_idxs)

             def compute_offset(ref_seq):
                 offset=0
                 for c in ref_seq:
                     if c in set(['-']):
                         offset+=1
                     else:
                         return offset

             #compute positions relative to alignmnet
             df_needle_alignment['ref_positions']=df_needle_alignment['ref_seq'].apply(compute_ref_positions)
             df_needle_alignment['read_offset']=df_needle_alignment['ref_seq'].apply(compute_offset)


             #INITIALIZATIONS
             re_find_indels=re.compile("(-*-)")
             re_find_indels_start=re.compile("(^-*-)")
             re_find_indels_end=re.compile("(-*-$)")
             re_find_substitutions=re.compile("(\.*\.)")

             effect_vector_insertion=np.zeros(len_amplicon)
             effect_vector_deletion=np.zeros(len_amplicon)
             effect_vector_mutation=np.zeros(len_amplicon)
             effect_vector_any=np.zeros(len_amplicon)

             effect_vector_insertion_mixed=np.zeros(len_amplicon)
             effect_vector_deletion_mixed=np.zeros(len_amplicon)
             effect_vector_mutation_mixed=np.zeros(len_amplicon)

             effect_vector_insertion_hdr=np.zeros(len_amplicon)
             effect_vector_deletion_hdr=np.zeros(len_amplicon)
             effect_vector_mutation_hdr=np.zeros(len_amplicon)

             effect_vector_insertion_noncoding=np.zeros(len_amplicon)
             effect_vector_deletion_noncoding=np.zeros(len_amplicon)
             effect_vector_mutation_noncoding=np.zeros(len_amplicon)

             hist_inframe=defaultdict(lambda :0)
             hist_frameshift=defaultdict(lambda :0)

             avg_vector_del_all=np.zeros(len_amplicon)
             avg_vector_ins_all=np.zeros(len_amplicon)

             print('SETUP for quant...')
             #look around the sgRNA(s) only?
             if cut_points and args.window_around_sgrna>0:
                include_idxs=[]
                half_window=max(1,args.window_around_sgrna/2)
                for cut_p in cut_points:
                    st=max(0,cut_p-half_window+1)
                    en=min(len(args.amplicon_seq)-1,cut_p+half_window+1)
                    include_idxs.append(range(st,en))
             else:
                include_idxs=range(len(args.amplicon_seq))

             print('about to do exclude idxs for quant...')
             exclude_idxs=[]

             if args.exclude_bp_from_left:
                exclude_idxs+=range(args.exclude_bp_from_left)

             if args.exclude_bp_from_right:
                exclude_idxs+=range(len_amplicon)[-args.exclude_bp_from_right:]


             #flatten the arrays to avoid errors with old numpy library
             include_idxs=np.ravel(include_idxs)
             exclude_idxs=np.ravel(exclude_idxs)

             include_idxs=set(np.setdiff1d(include_idxs,exclude_idxs))



             #handy generator to split in chunks the dataframe, np.split_array is slow!
             def get_chunk(df_needle_alignment,n_processes=args.n_processes):
                 for g,df in df_needle_alignment.groupby(np.arange(len(df_needle_alignment)) // (len(df_needle_alignment)/(args.n_processes-1))):
                     yield df



             #Use a Pool of processes, or just a single process
             if args.n_processes > 1:
                info('[CRISPResso quantification is running in parallel mode with %d processes]' % min(df_needle_alignment.shape[0],args.n_processes) )
                pool = mp.Pool(processes=min(df_needle_alignment.shape[0],args.n_processes))
                chunks_computed=[]
                for result in pool.imap(process_df_chunk,get_chunk(df_needle_alignment)):
                     df_needle_alignment_chunk, effect_vector_insertion_chunk,effect_vector_deletion_chunk,\
                     effect_vector_mutation_chunk,effect_vector_any_chunk,effect_vector_insertion_mixed_chunk,effect_vector_deletion_mixed_chunk,\
                     effect_vector_mutation_mixed_chunk,effect_vector_insertion_hdr_chunk,effect_vector_deletion_hdr_chunk,effect_vector_mutation_hdr_chunk,\
                     effect_vector_insertion_noncoding_chunk,effect_vector_deletion_noncoding_chunk,effect_vector_mutation_noncoding_chunk,hist_inframe_chunk,\
                     hist_frameshift_chunk,avg_vector_del_all_chunk,avg_vector_ins_all_chunk,MODIFIED_FRAMESHIFT_chunk,MODIFIED_NON_FRAMESHIFT_chunk,NON_INDELS_NON_FRAMESHIFT_chunk,\
                     SPLICING_SITES_MODIFIED_chunk=result

                     chunks_computed.append(df_needle_alignment_chunk)
                     effect_vector_insertion+=effect_vector_insertion_chunk
                     effect_vector_deletion+=effect_vector_deletion_chunk
                     effect_vector_mutation+=effect_vector_mutation_chunk
                     effect_vector_any+=effect_vector_any_chunk
                     effect_vector_insertion_mixed+=effect_vector_insertion_mixed_chunk
                     effect_vector_deletion_mixed+=effect_vector_deletion_mixed_chunk
                     effect_vector_mutation_mixed+=effect_vector_mutation_mixed_chunk
                     effect_vector_insertion_hdr+=effect_vector_insertion_hdr_chunk
                     effect_vector_deletion_hdr+=effect_vector_deletion_hdr_chunk
                     effect_vector_mutation_hdr+=effect_vector_mutation_hdr_chunk
                     effect_vector_insertion_noncoding+=effect_vector_insertion_noncoding_chunk
                     effect_vector_deletion_noncoding+=effect_vector_deletion_noncoding_chunk
                     effect_vector_mutation_noncoding+=effect_vector_mutation_noncoding_chunk
                     add_hist(hist_inframe_chunk,hist_inframe)
                     add_hist(hist_frameshift_chunk,hist_frameshift)
                     avg_vector_del_all+=avg_vector_del_all_chunk
                     avg_vector_ins_all+=avg_vector_ins_all_chunk
                     MODIFIED_FRAMESHIFT+=MODIFIED_FRAMESHIFT_chunk
                     MODIFIED_NON_FRAMESHIFT+=MODIFIED_NON_FRAMESHIFT_chunk
                     NON_INDELS_NON_FRAMESHIFT+=NON_INDELS_NON_FRAMESHIFT_chunk
                     SPLICING_SITES_MODIFIED+=SPLICING_SITES_MODIFIED_chunk

                pool.close()
                pool.join()
                df_needle_alignment=pd.concat(chunks_computed)
                del chunks_computed

             else:
                 df_needle_alignment, effect_vector_insertion,\
                 effect_vector_deletion,effect_vector_mutation,\
                 effect_vector_any,effect_vector_insertion_mixed,\
                 effect_vector_deletion_mixed,effect_vector_mutation_mixed,\
                 effect_vector_insertion_hdr,effect_vector_deletion_hdr,\
                 effect_vector_mutation_hdr,effect_vector_insertion_noncoding,\
                 effect_vector_deletion_noncoding,effect_vector_mutation_noncoding,\
                 hist_inframe,hist_frameshift,avg_vector_del_all,\
                 avg_vector_ins_all,MODIFIED_FRAMESHIFT,MODIFIED_NON_FRAMESHIFT,\
                 NON_INDELS_NON_FRAMESHIFT,SPLICING_SITES_MODIFIED= process_df_chunk(df_needle_alignment)
             print('finished process_df_chunk...')


             N_INDELS=df_needle_alignment['NHEJ'].sum()
             N_UNMODIFIED=df_needle_alignment['UNMODIFIED'].sum()
             #N_PARTIAL=df_needle_alignment['PARTIAL'].sum()
             N_MIXED_HDR_NHEJ=df_needle_alignment['MIXED'].sum()
             N_REPAIRED=df_needle_alignment['HDR'].sum()
             N_CUTSUBS=df_needle_alignment['CUTSUB'].sum()
             #import pdb; pdb.set_trace()

             #disable known division warning
             with np.errstate(divide='ignore',invalid='ignore'):

                 effect_vector_combined=100*effect_vector_any/float(N_ALIGNED)

                 avg_vector_ins_all/=(effect_vector_insertion+effect_vector_insertion_hdr+effect_vector_insertion_mixed)
                 avg_vector_del_all/=(effect_vector_deletion+effect_vector_deletion_hdr+effect_vector_deletion_mixed)

             avg_vector_ins_all[np.isnan(avg_vector_ins_all)]=0
             avg_vector_del_all[np.isnan(avg_vector_del_all)]=0
             avg_vector_ins_all[np.isinf(avg_vector_ins_all)]=0
             avg_vector_del_all[np.isinf(avg_vector_del_all)]=0

             if PERFORM_FRAMESHIFT_ANALYSIS:
                 if not dict(hist_inframe):
                    hist_inframe={0:0}

                 if not dict(hist_frameshift):
                    hist_frameshift={0:0}

             print('==>Done.')

             print('Calculating indel distribution based on the length of the reads...')

             df_needle_alignment['effective_len']=df_needle_alignment.apply(lambda row:  len_amplicon+row.n_inserted-row.n_deleted,axis=1)

             print('==>Done.')

             #write alleles table
             print ('Calculating alleles frequencies...')

             def get_ref_positions(row,df_alignment):
                return list(df_alignment.ix[(row.Aligned_Sequence ,row.Reference_Sequence),'ref_positions'][0])

             df_alleles=df_needle_alignment.groupby(['align_seq','ref_seq','NHEJ','UNMODIFIED','HDR','PARTIAL','CUTSUB','n_deleted','n_inserted','n_mutated',]).size()
             #df_alleles=df_needle_alignment.groupby(['align_seq','ref_seq','NHEJ','UNMODIFIED','HDR','SUB','n_deleted','n_inserted','n_mutated',]).size()
             df_alleles=df_alleles.reset_index()
             df_alleles.rename(columns={0:'#Reads','align_seq':'Aligned_Sequence','ref_seq':'Reference_Sequence'},inplace=True)
             #df_alleles.set_index('Aligned_Sequence',inplace=True)
             df_alleles['%Reads']=df_alleles['#Reads']/df_alleles['#Reads'].sum()*100

             if np.sum(np.array(map(int,pd.__version__.split('.')))*(100,10,1))< 170:
                df_alleles.sort('#Reads',ascending=False,inplace=True)
             else:
                df_alleles.sort_values(by='#Reads',ascending=False,inplace=True)

             
             #add ref positions for the plot around the cut sites
             df_needle_alignment.set_index(['align_seq','ref_seq'],inplace=True)
             df_needle_alignment.sort_index(inplace=True)
	     #import pdb; pdb.set_trace()
             #df_alleles.head(100).to_csv('df_alleles.csv')
             #df_needle_alignment.head(100).to_csv('df_needle_alignment.csv')
             df_alleles['ref_positions']=df_alleles.apply(lambda x: get_ref_positions(x,df_needle_alignment),axis=1).values

             print('==>Done.')


             print('Making Plots...')
             #plot effective length
             if args.guide_seq:
                 min_cut=min(cut_points)
                 max_cut=max(cut_points)
                 xmin,xmax=-min_cut,len_amplicon-max_cut
             else:
                 min_cut=len_amplicon/2
                 max_cut=len_amplicon/2
                 xmin,xmax=-min_cut,+max_cut


             hdensity,hlengths=np.histogram(df_needle_alignment.effective_len-len_amplicon,np.arange(xmin,xmax))
             hlengths=hlengths[:-1]
             center_index=np.nonzero(hlengths==0)[0][0]

             fig=plt.figure(figsize=(8.3,8))

             plt.bar(0,hdensity[center_index],color='red',linewidth=0)
             #plt.hold(True)
             barlist=plt.bar(hlengths,hdensity,align='center',linewidth=0)
             barlist[center_index].set_color('r')
             plt.xlim([xmin,xmax])
             plt.ylabel('Sequences (no.)')
             plt.xlabel('Indel size (bp)')
             plt.ylim([0,hdensity.max()*1.2])
             plt.title('Indel size distribution')
             lgd=plt.legend(['No indel','Indel'],loc='center', bbox_to_anchor=(0.5, -0.22),ncol=1, fancybox=True, shadow=True)
             #lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.28),ncol=1, fancybox=True, shadow=True)
             lgd.legendHandles[0].set_height(3)
             lgd.legendHandles[1].set_height(3)
             plt.savefig(_jp('1a.Indel_size_distribution_n_sequences.pdf'),bbox_inches='tight')
             if args.save_also_png:
                     plt.savefig(_jp('1a.Indel_size_distribution_n_sequences.png'),bbox_inches='tight')


             plt.figure(figsize=(8.3,8))
             plt.bar(0,hdensity[center_index]/(float(hdensity.sum()))*100.0,color='red',linewidth=0)
             #plt.hold(True)
             barlist=plt.bar(hlengths,hdensity/(float(hdensity.sum()))*100.0,align='center',linewidth=0)
             barlist[center_index].set_color('r')
             plt.xlim([xmin,xmax])
             plt.title('Indel size distribution')
             plt.ylabel('Sequences (%)')
             plt.xlabel('Indel size (bp)')
             #lgd=plt.legend(['No indel','Indel'])
             lgd=plt.legend(['No indel','Indel'],loc='center', bbox_to_anchor=(0.5, -0.22),ncol=1, fancybox=True, shadow=True)
             lgd.legendHandles[0].set_height(3)
             lgd.legendHandles[1].set_height(3)
             print('Plot 1a done...')

             plt.savefig(_jp('1b.Indel_size_distribution_percentage.pdf'),bbox_inches='tight')
             if args.save_also_png:
                     plt.savefig(_jp('1b.Indel_size_distribution_percentage.png'),bbox_inches='tight')

             print('Plot 1b done...')

             ####PIE CHARTS FOR HDR/NHEJ/MIXED/EVENTS###

             if args.expected_hdr_amplicon_seq:
                 fig=plt.figure(figsize=(12*1.5,14.5*1.5))
                 ax1 = plt.subplot2grid((6,3), (0, 0), colspan=3, rowspan=5)
                 patches, texts, autotexts =ax1.pie([N_UNMODIFIED,N_CUTSUBS,N_INDELS,N_REPAIRED],\
                                                                   labels=['Unmodified\n(%d reads)' %N_UNMODIFIED,\
                                                                           'Substitutions\n(%d reads)' %N_CUTSUBS,
                                                                           'NHEJ\n(%d reads)' % N_INDELS, \
                                                                           'HDR\n(%d reads)' %N_REPAIRED,
                                                                           ],\
                                                                   explode=(0,0,0,0),\
                                                                   colors=[(1,0,0,0.2),(0,1,1,0.2),(0,0,1,0.2),(0,1,0,0.2)],autopct='%1.1f%%')

                 if cut_points or args.donor_seq:
                    ax2 = plt.subplot2grid((6,3), (5, 0), colspan=3, rowspan=1)
                    ax2.plot([0,len_amplicon],[0,0],'-k',lw=2,label='Amplicon sequence')
                    #plt.hold(True)

                    if args.donor_seq:
                        ax2.plot(core_donor_seq_st_en,[0,0],'-',lw=10,c=(0,1,0,0.5),label='Donor Sequence')

                    if cut_points:
                        ax2.plot(cut_points+offset_plots,np.zeros(len(cut_points)),'vr', ms=24,label='Predicted Cas9 cleavage site/s')

                    for idx,sgRNA_int in enumerate(sgRNA_intervals):
                         if idx==0:
                            ax2.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA')
                         else:
                            ax2.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_')

                    plt.legend(bbox_to_anchor=(0, 0, 1., 0),  ncol=1, mode="expand", borderaxespad=0.,numpoints=1)
                    plt.xlim(0,len_amplicon)
                    plt.axis('off')

                 proptease = fm.FontProperties()
                 proptease.set_size('xx-large')
                 plt.setp(autotexts, fontproperties=proptease)
                 plt.setp(texts, fontproperties=proptease)
                 plt.savefig(_jp('2.Unmodified_NHEJ_HDR_pie_chart.pdf'),pad_inches=1,bbox_inches='tight')
                 if args.save_also_png:
                         plt.savefig(_jp('2.Unmodified_NHEJ_HDR_pie_chart.png'),pad_inches=1,bbox_inches='tight')
                 print('Plot 2 done...')


             else:
                 fig=plt.figure(figsize=(12*1.5,14.5*1.5))
                 ax1 = plt.subplot2grid((6,3), (0, 0), colspan=3, rowspan=5)
                 patches, texts, autotexts =ax1.pie([N_UNMODIFIED/N_ALIGNED*100,N_INDELS/N_ALIGNED*100],\
                                                   labels=['Unmodified\n(%d reads)' %N_UNMODIFIED,\
                                                           'NHEJ\n(%d reads)' % N_INDELS],\
                                                   explode=(0,0),colors=[(1,0,0,0.2),(0,0,1,0.2)],autopct='%1.1f%%')

                 if cut_points:
                    ax2 = plt.subplot2grid((6,3), (5, 0), colspan=3, rowspan=1)
                    ax2.plot([0,len_amplicon],[0,0],'-k',lw=2,label='Amplicon sequence')
                    #plt.hold(True)


                    for idx,sgRNA_int in enumerate(sgRNA_intervals):
                         if idx==0:
                            ax2.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                         else:
                            ax2.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')

                    ax2.plot(cut_points+offset_plots,np.zeros(len(cut_points)),'vr', ms=12,label='Predicted Cas9 cleavage site/s')
                    plt.legend(bbox_to_anchor=(0, 0, 1., 0),  ncol=1, mode="expand", borderaxespad=0.,numpoints=1,prop={'size':'large'})
                    plt.xlim(0,len_amplicon)
                    plt.axis('off')

                 proptease = fm.FontProperties()
                 proptease.set_size('xx-large')
                 plt.setp(autotexts, fontproperties=proptease)
                 plt.setp(texts, fontproperties=proptease)
                 plt.savefig(_jp('2.Unmodified_NHEJ_pie_chart.pdf'),pad_inches=1,bbox_inches='tight')
                 if args.save_also_png:
                         plt.savefig(_jp('2.Unmodified_NHEJ_pie_chart.png'),pad_inches=1,bbox_inches='tight')
                 print('Plot 2 done...')


             ###############################################################################################################################################


             #(3) a graph of frequency of deletions and insertions of various sizes (deletions could be consider as negative numbers and insertions as positive);


             def calculate_range(df,column_name):
                df_not_zero=df.ix[df[column_name]>0,column_name]
                try:
                    r=max(15,int(np.round(np.percentile(df_not_zero,99))))
                except:
                    r=15
                return r

             range_mut=calculate_range(df_needle_alignment,'n_mutated')
             range_ins=calculate_range(df_needle_alignment,'n_inserted')
             range_del=calculate_range(df_needle_alignment,'n_deleted')

             y_values_mut,x_bins_mut=plt.histogram(df_needle_alignment['n_mutated'],bins=range(0,range_mut))
             y_values_ins,x_bins_ins=plt.histogram(df_needle_alignment['n_inserted'],bins=range(0,range_ins))
             y_values_del,x_bins_del=plt.histogram(df_needle_alignment['n_deleted'],bins=range(0,range_del))

             fig=plt.figure(figsize=(26,6.5))


             ax=fig.add_subplot(1,3,1)
             ax.bar(x_bins_ins[:-1],y_values_ins,align='center',linewidth=0,color=(0,0,1))
             barlist=ax.bar(x_bins_ins[:-1],y_values_ins,align='center',linewidth=0,color=(0,0,1))
             barlist[0].set_color('r')

             plt.title('Insertions')
             plt.xlabel('Size (bp)')
             plt.ylabel('Sequences % (no.)')
             lgd=plt.legend(['Non-insertion','Insertion'][::-1], bbox_to_anchor=(.82, -0.22),ncol=1, fancybox=True, shadow=True)
             lgd.legendHandles[0].set_height(6)
             lgd.legendHandles[1].set_height(6)
             plt.xlim(xmin=-1)
             y_label_values= np.round(np.linspace(0, min(N_ALIGNED,max(ax.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
             plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/N_ALIGNED*100,n_reads) for n_reads in y_label_values])

             ax=fig.add_subplot(1,3,2)
             ax.bar(-x_bins_del[:-1],y_values_del,align='center',linewidth=0,color=(0,0,1))
             barlist=ax.bar(-x_bins_del[:-1],y_values_del,align='center',linewidth=0,color=(0,0,1))
             barlist[0].set_color('r')
             plt.title('Deletions')
             plt.xlabel('Size (bp)')
             plt.ylabel('Sequences % (no.)')
             lgd=plt.legend(['Non-deletion','Deletion'][::-1], bbox_to_anchor=(.82, -0.22),ncol=1, fancybox=True, shadow=True)
             lgd.legendHandles[0].set_height(6)
             lgd.legendHandles[1].set_height(6)
             plt.xlim(xmax=1)
             y_label_values= np.round(np.linspace(0, min(N_ALIGNED,max(ax.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
             plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/N_ALIGNED*100,n_reads) for n_reads in y_label_values])

             ax=fig.add_subplot(1,3,3)
             ax.bar(x_bins_mut[:-1],y_values_mut,align='center',linewidth=0,color=(0,0,1))
             barlist=ax.bar(x_bins_mut[:-1],y_values_mut,align='center',linewidth=0,color=(0,0,1))
             barlist[0].set_color('r')
             plt.title('Substitutions')
             plt.xlabel('Positions substituted (number)')
             plt.ylabel('Sequences % (no.)')
             lgd=plt.legend(['Non-substitution','Substitution'][::-1] ,bbox_to_anchor=(.82, -0.22),ncol=1, fancybox=True, shadow=True)
             lgd.legendHandles[0].set_height(6)
             lgd.legendHandles[1].set_height(6)
             plt.xlim(xmin=-1)
             y_label_values= np.round(np.linspace(0, min(N_ALIGNED,max(ax.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
             plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/N_ALIGNED*100,n_reads) for n_reads in y_label_values])
             plt.tight_layout()

             plt.savefig(_jp('3.Insertion_Deletion_Substitutions_size_hist.pdf'),bbox_inches='tight')
             if args.save_also_png:
                     plt.savefig(_jp('3.Insertion_Deletion_Substitutions_size_hist.png'),bbox_inches='tight')
             print('Plot 3 done...')

             #(4) another graph with the frequency that each nucleotide within the amplicon was modified in any way (perhaps would consider insertion as modification of the flanking nucleotides);

             #Indels location Plots

             plt.figure(figsize=(10,10))

             y_max=max(effect_vector_any)*1.2

             plt.plot(effect_vector_any,'r',lw=3,label='Combined Insertions/Deletions/Substitutions')
             #plt.hold(True)

             if cut_points:

                 for idx,cut_point in enumerate(cut_points):
                     if idx==0:
                             plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                     else:
                             plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='_nolegend_')


                 for idx,sgRNA_int in enumerate(sgRNA_intervals):
                     if idx==0:
                        plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                     else:
                        plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')


             lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.23),ncol=1, fancybox=True, shadow=True)
             y_label_values= np.round(np.linspace(0, min(N_ALIGNED,max(ax.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
             plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/float(N_ALIGNED)*100, n_reads) for n_reads in y_label_values])
             plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )

             plt.title('Mutation position distribution')
             plt.xlabel('Reference amplicon position (bp)')
             plt.ylabel('Sequences % (no.)')
             plt.ylim(0,max(1,y_max))
             print('Fig 4 done...')
             plt.xlim(xmax=len(args.amplicon_seq)-1)
             plt.savefig(_jp('4a.Combined_Insertion_Deletion_Substitution_Locations.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
             if args.save_also_png:
                     plt.savefig(_jp('4a.Combined_Insertion_Deletion_Substitution_Locations.png'),bbox_extra_artists=(lgd,), bbox_inches='tight',pad=1)
             print('Plot 4a done...')

#-----------------------------------------------------------------------------------------------------------
             #NHEJ
             plt.figure(figsize=(10,10))
             #effect_vector_insertion[0] = 0
             #effect_vector_insertion[300] = 0
             #plt.plot(effect_vector_insertion[0:300],'r',lw=3,label='Insertions')
             plt.plot(effect_vector_insertion,'r',lw=3,label='Insertions')
             #plt.hold(True)
             plt.plot(effect_vector_deletion,'m',lw=3,label='Deletions')
             plt.plot(effect_vector_mutation,'g',lw=3,label='Substitutions')

             y_max=max(max(effect_vector_insertion),max(effect_vector_deletion),max(effect_vector_mutation))*1.2


             if cut_points:

                 for idx,cut_point in enumerate(cut_points):
                     if idx==0:
                             plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                     else:
                             plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='_nolegend_')


                 #for idx,sgRNA_int in enumerate(sgRNA_intervals):
                 #    if idx==0:
                 #       plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                 #    else:
                 #       plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')

             lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.28),ncol=1, fancybox=True, shadow=True)
             y_label_values=np.arange(0,y_max,max(y_max/6.0,1)) #SKW fix div 0 rounding error
             plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_ALIGNED)*100,n_reads/float(N_INDELS)*100, n_reads) for n_reads in y_label_values])
             plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )

             plt.xlabel('Reference amplicon position (bp)')
             plt.ylabel('Sequences: % Total ( % NHEJ, no. )')
             plt.ylim(0,max(1,y_max))
             plt.xlim(xmax=len(args.amplicon_seq)-1)

             plt.title('Mutation position distribution of NHEJ')
             plt.savefig(_jp('4b.NHEJ_%s.pdf' % args.name),bbox_extra_artists=(lgd,), bbox_inches='tight')
             if args.save_also_png:
                     plt.savefig(_jp('4b.Insertion_Deletion_Substitution_Locations_NHEJ.png'),bbox_extra_artists=(lgd,), bbox_inches='tight',pad=1)
             print('Plot 4b done...')


#-----------------------------------------------------------------------------------------------------------


             if args.expected_hdr_amplicon_seq:

                 #HDR
                 plt.figure(figsize=(10,10))
                 plt.plot(effect_vector_insertion_hdr,'r',lw=3,label='Insertions')
                 #plt.hold(True)
                 plt.plot(effect_vector_deletion_hdr,'m',lw=3,label='Deletions')
                 plt.plot(effect_vector_mutation_hdr,'g',lw=3,label='Substitutions')

                 y_max=max(max(effect_vector_insertion_hdr),max(effect_vector_deletion_hdr),max(effect_vector_mutation_hdr))*1.2

                 if cut_points:

                         for idx,cut_point in enumerate(cut_points):
                             if idx==0:
                                     plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                             else:
                                     plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='_nolegend_')


                         for idx,sgRNA_int in enumerate(sgRNA_intervals):
                             if idx==0:
                                plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                             else:
                                plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')


                 lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.28),ncol=1, fancybox=True, shadow=True)
                 y_label_values=np.arange(0,y_max,max(y_max/6,1)).astype(int) #SKW fix div 0 rounding error
                 plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_ALIGNED)*100,n_reads/float(N_REPAIRED)*100, n_reads) for n_reads in y_label_values])
                 plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )

                 plt.xlabel('Reference amplicon position (bp)')
                 plt.ylabel('Sequences: % Total ( % HDR, no. )')
                 plt.ylim(0,max(1,y_max))
                 plt.xlim(xmax=len(args.amplicon_seq)-1)
                 plt.title('Mutation position distribution of HDR')
                 plt.savefig(_jp('4c.Insertion_Deletion_Substitution_Locations_HDR.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
                 if args.save_also_png:
                     plt.savefig(_jp('4c.Insertion_Deletion_Substitution_Locations_HDR.png'),bbox_extra_artists=(lgd,), bbox_inches='tight',pad=1)
                 print('Plot 4c done...')

#-----------------------------------------------------------------------------------------------------------

                 #MIXED
                 plt.figure(figsize=(10,10))
                 plt.plot(effect_vector_insertion_mixed,'r',lw=3,label='Insertions')
                 #plt.hold(True)
                 plt.plot(effect_vector_deletion_mixed,'m',lw=3,label='Deletions')
                 plt.plot(effect_vector_mutation_mixed,'g',lw=3,label='Substitutions')

                 y_max=max(max(effect_vector_insertion_mixed),max(effect_vector_deletion_mixed),max(effect_vector_mutation_mixed))*1.2

                 if cut_points:

                         for idx,cut_point in enumerate(cut_points):
                             if idx==0:
                                     plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                             else:
                                     plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='_nolegend_')

                         for idx,sgRNA_int in enumerate(sgRNA_intervals):
                             if idx==0:
                                plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                             else:
                                plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')

                 lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.28),ncol=1, fancybox=True, shadow=True)
                 y_label_values=np.arange(0,y_max,max(y_max/6.0,1))
                 plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_ALIGNED)*100,n_reads/float(N_MIXED_HDR_NHEJ)*100, n_reads) for n_reads in y_label_values])
                 # plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )
                 plt.xticks(np.linspace(1,len_amplicon,6).astype(int) )  # SKW from arange->linspace

                 plt.xlabel('Reference amplicon position (bp)')
                 plt.ylabel('Sequences: % Total ( % mixed HDR-NHEJ, no. )')
                 plt.ylim(0,max(1,y_max))
                 plt.xlim(xmax=len(args.amplicon_seq)-1)
                 plt.title('Mutation position distribution of mixed HDR-NHEJ')
                 plt.savefig(_jp('4d.Insertion_Deletion_Substitution_Locations_Mixed_HDR_NHEJ.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
                 if args.save_also_png:
                         plt.savefig(_jp('4d.Insertion_Deletion_Substitution_Locations_Mixed_HDR_NHEJ.png'),bbox_extra_artists=(lgd,), bbox_inches='tight',pad=1)
                 print('Plot 4d done...')
#-----------------------------------------------------------------------------------------------------------


            #Position dependent indels plot
             fig=plt.figure(figsize=(24,10))
             ax1=fig.add_subplot(1,2,1)
             markerline, stemlines, baseline=ax1.stem(avg_vector_ins_all,'r',lw=3,markerfmt="s",markerline=None,s=50)
             plt.setp(markerline, 'markerfacecolor', 'r', 'markersize', 8)
             plt.setp(baseline, 'linewidth', 0)
             plt.setp(stemlines, 'color', 'r','linewidth',3)
             #plt.hold(True)
             y_max=max(avg_vector_ins_all)*1.2
             if cut_points:

                 for idx,cut_point in enumerate(cut_points):
                     if idx==0:
                             ax1.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                     else:
                             ax1.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='_nolegend_')

             #plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )
             plt.xticks(np.linspace(1,len_amplicon,6).astype(int) )  # SKW from arange->linspace
             plt.xlabel('Reference amplicon position (bp)')
             plt.ylabel('Average insertion length')
             plt.ylim(0,max(1,y_max))
             plt.xlim(xmax=len_amplicon-1)
             ax1.set_title('Position dependent insertion size')
             plt.tight_layout()

             ax2=fig.add_subplot(1,2,2)
             markerline, stemlines, baseline=ax2.stem(avg_vector_del_all,'r',lw=3,markerfmt="s",markerline=None,s=50)
             plt.setp(markerline, 'markerfacecolor', 'm', 'markersize', 8)
             plt.setp(baseline, 'linewidth', 0)
             plt.setp(stemlines, 'color', 'm','linewidth',3)
             #plt.hold(True)
             y_max=max(avg_vector_del_all)*1.2
             if cut_points:

                 for idx,cut_point in enumerate(cut_points):
                     if idx==0:
                             ax2.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                     else:
                             ax2.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='_nolegend_')

             #plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )
             plt.xticks(np.linspace(1,len_amplicon,6).astype(int) )  # SKW from arange->linspace
             plt.xlabel('Reference amplicon position (bp)')
             plt.ylabel('Average deletion length')

             plt.ylim(ymin=0,ymax=max(1,y_max))
             plt.xlim(xmax=len_amplicon-1)
             ax2.set_title('Position dependent deletion size')

             plt.tight_layout()


             plt.savefig(_jp('4e.Position_dependent_average_indel_size.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
             if args.save_also_png:
                 plt.savefig(_jp('4e.Position_dependent_average_indel_size.png'),bbox_extra_artists=(lgd,), bbox_inches='tight')
             print('Plot 4e done...')
#-----------------------------------------------------------------------------------------------------------


             if PERFORM_FRAMESHIFT_ANALYSIS:
                 #make frameshift plots
                 fig=plt.figure(figsize=(12*1.5,14.5*1.5))
                 ax1 = plt.subplot2grid((6,3), (0, 0), colspan=3, rowspan=5)
                 patches, texts, autotexts =ax1.pie([MODIFIED_FRAMESHIFT,\
                                                    MODIFIED_NON_FRAMESHIFT,\
                                                    NON_INDELS_NON_FRAMESHIFT],\
                                                    labels=['Frameshift mutation\n(%d reads)' %MODIFIED_FRAMESHIFT,\
                                                           'In-frame mutation\n(%d reads)' % MODIFIED_NON_FRAMESHIFT,\
                                                           'Noncoding mutation\n(%d reads)' %NON_INDELS_NON_FRAMESHIFT],\
                                                    explode=(0.0,0.0,0.0),\
                                                    colors=[(0.89019608,  0.29019608,  0.2, 0.8),(0.99215686,  0.73333333,  0.51764706,0.8),(0.99607843,  0.90980392,  0.78431373,0.8)],\
                                                    autopct='%1.1f%%')

                 ax2 = plt.subplot2grid((6,3), (5, 0), colspan=3, rowspan=1)
                 ax2.plot([0,len_amplicon],[0,0],'-k',lw=2,label='Amplicon sequence')
                 #plt.hold(True)

                 for idx,exon_interval in enumerate(exon_intervals):
                     if idx==0:
                         ax2.plot(exon_interval,[0,0],'-',lw=10,c=(0,0,1,0.5),label='Coding sequence/s',solid_capstyle='butt')
                     else:
                         ax2.plot(exon_interval,[0,0],'-',lw=10,c=(0,0,1,0.5),label='_nolegend_',solid_capstyle='butt')

                 if cut_points:
                    ax2.plot(cut_points+offset_plots,np.zeros(len(cut_points)),'vr', ms=25,label='Predicted Cas9 cleavage site/s')

                 plt.legend(bbox_to_anchor=(0, 0, 1., 0),  ncol=1, mode="expand", borderaxespad=0.,numpoints=1)
                 plt.xlim(0,len_amplicon)
                 plt.axis('off')

                 proptease = fm.FontProperties()
                 proptease.set_size('xx-large')
                 plt.setp(autotexts, fontproperties=proptease)
                 plt.setp(texts, fontproperties=proptease)
                 plt.savefig(_jp('5.Frameshift_In-frame_mutations_pie_chart.pdf'),pad_inches=1,bbox_inches='tight')
                 if args.save_also_png:
                         plt.savefig(_jp('5.Frameshift_In-frame_mutations_pie_chart.png'),pad_inches=1,bbox_inches='tight')
                 print('Plot 5 done...')


                 #profiles-----------------------------------------------------------------------------------
                 fig=plt.figure(figsize=(22,10))
                 ax1=fig.add_subplot(2,1,1)
                 x,y=map(np.array,zip(*[a for a in hist_frameshift.iteritems()]))
                 y=y/float(sum(hist_frameshift.values()))*100
                 ax1.bar(x-0.5,y)
                 ax1.set_xlim(-30.5,30.5)
                 ax1.set_frame_on(False)
                 ax1.set_xticks([idx for idx in range(-30,31) if idx % 3])
                 ax1.tick_params(which='both',      # both major and minor ticks are affected
                    bottom='off',      # ticks along the bottom edge are off
                    top='off',         # ticks along the top edge are off
                    labelbottom='on') # labels along the bottom edge are off)
                 ax1.yaxis.tick_left()
                 xmin, xmax = ax1.get_xaxis().get_view_interval()
                 ymin, ymax = ax1.get_yaxis().get_view_interval()
                 ax1.set_xticklabels([str(idx)  for idx in [idx for idx in range(-30,31) if idx % 3]],rotation='vertical')
                 plt.title('Frameshift profile')
                 ax1.tick_params(axis='both', which='major', labelsize=32)
                 ax1.tick_params(axis='both', which='minor', labelsize=32)
                 plt.tight_layout()
                 plt.ylabel('%')

                 ax2=fig.add_subplot(2,1,2)
                 x,y=map(np.array,zip(*[a for a in hist_inframe.iteritems()]))
                 y=y/float(sum(hist_inframe.values()))*100
                 ax2.bar(x-0.5,y,color=(0,1,1,0.2))
                 ax2.set_xlim(-30.5,30.5)
                 ax2.set_frame_on(False)
                 ax2.set_xticks([idx for idx in range(-30,31) if (idx % 3 ==0) ])
                 ax2.tick_params(which='both',      # both major and minor ticks are affected
                    bottom='off',      # ticks along the bottom edge are off
                    top='off',         # ticks along the top edge are off
                    labelbottom='on') # labels along the bottom edge are off)
                 ax2.yaxis.tick_left()
                 xmin, xmax = ax2.xaxis.get_view_interval()
                 ymin, ymax = ax2.yaxis.get_view_interval()
                 ax2.set_xticklabels([str(idx)  for idx in [idx for idx in range(-30,31) if (idx % 3==0)]],rotation='vertical')
                 plt.title('In-frame profile')
                 plt.tight_layout()
                 plt.ylabel('%')
                 ax2.tick_params(axis='both', which='major', labelsize=32)
                 ax2.tick_params(axis='both', which='minor', labelsize=32)
                 plt.tight_layout()

                 plt.savefig(_jp('6.Frameshift_In-frame_mutation_profiles.pdf'),pad_inches=1,bbox_inches='tight')
                 if args.save_also_png:
                     plt.savefig(_jp('6.Frameshift_In-frame_mutation_profiles.png'),pad_inches=1,bbox_inches='tight')
                 print('Plot 6 done...')

                 #-----------------------------------------------------------------------------------------------------------
                 fig=plt.figure(figsize=(12*1.5,12*1.5))
                 ax=fig.add_subplot(1,1,1)
                 patches, texts, autotexts =ax.pie([SPLICING_SITES_MODIFIED,\
                                                   (df_needle_alignment.shape[0] - SPLICING_SITES_MODIFIED)],\
                                                   labels=['Potential splice sites modified\n(%d reads)' %SPLICING_SITES_MODIFIED,\
                                                           'Unmodified\n(%d reads)' % (df_needle_alignment.shape[0]- SPLICING_SITES_MODIFIED)],\
                                                   explode=(0.0,0),\
                                                   colors=[(0.89019608,  0.29019608,  0.2, 0.8),(0.99607843,  0.90980392,  0.78431373,0.8)],\
                                                   autopct='%1.1f%%')
                 proptease = fm.FontProperties()
                 proptease.set_size('xx-large')
                 plt.setp(autotexts, fontproperties=proptease)
                 plt.setp(texts, fontproperties=proptease)
                 plt.savefig(_jp('8.Potential_Splice_Sites_pie_chart.pdf'),pad_inches=1,bbox_inches='tight')
                 if args.save_also_png:
                     plt.savefig(_jp('8.Potential_Splice_Sites_pie_chart.png'),pad_inches=1,bbox_inches='tight')
                 print('Plot 8 done...')

                 #non coding
                 plt.figure(figsize=(10,10))
                 plt.plot(effect_vector_insertion_noncoding,'r',lw=3,label='Insertions')
                 #plt.hold(True)
                 plt.plot(effect_vector_deletion_noncoding,'m',lw=3,label='Deletions')
                 plt.plot(effect_vector_mutation_noncoding,'g',lw=3,label='Substitutions')

                 y_max=max(max(effect_vector_insertion_noncoding),max(effect_vector_deletion_noncoding),max(effect_vector_mutation_noncoding))*1.2


                 if cut_points:

                     for idx,cut_point in enumerate(cut_points):
                         if idx==0:
                                 plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                         else:
                                 plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='_nolegend_')

                         for idx,sgRNA_int in enumerate(sgRNA_intervals):
                             if idx==0:
                                plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                             else:
                                plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')

                 lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.28),ncol=1, fancybox=True, shadow=True)
                 #plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )
                 plt.xticks(np.linspace(1,len_amplicon,6).astype(int) )  # SKW from arange->linspace

                 plt.xlabel('Reference amplicon position (bp)')
                 plt.ylabel('Sequences (no.)')
                 plt.ylim(0,max(1,y_max))
                 plt.xlim(xmax=len(args.amplicon_seq)-1)
                 plt.title('Noncoding mutation position distribution')
                 plt.savefig(_jp('7.Insertion_Deletion_Substitution_Locations_Noncoding.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
                 if args.save_also_png:
                         plt.savefig(_jp('7.Insertion_Deletion_Substitution_Locations_Noncoding.png'),bbox_extra_artists=(lgd,), bbox_inches='tight')
                 print('Plot 7 done...')

             print('After FS analysis, before plots around cut sites...')
             ##new plots alleles around cut_sites

             for sgRNA,cut_point in zip(sgRNA_sequences,cut_points):
                 #print sgRNA,cut_point
                 len_amplicon = len(args.amplicon_seq)
                 #offset_to_plot = (min(len_amplicon - cut_point,cut_point)) - 20
                 distance_to_end = min(len_amplicon - cut_point,cut_point)
                 if distance_to_end < 75:
                      offset_to_plot = distance_to_end - 5
                 else:
                      offset_to_plot = 75
                 
                 # OFFSET WINDOW CHANGE HERE
                 # offset_to_plot = 37;
                 print('Printing from %d %d %d' % (offset_to_plot,len_amplicon,cut_point))
                 df_allele_around_cut=get_dataframe_around_cut(df_alleles, cut_point, offset_to_plot)
                 #df_allele_around_cut=get_dataframe_around_cut(df_alleles, cut_point,args.offset_around_cut_to_plot)

                 #write alleles table to file
                 df_allele_around_cut.to_csv(_jp('Alleles_frequency_table_around_cut_site_for_%s.txt' % sgRNA),sep='\t',header=True)

                 plot_alleles_table(offset_to_plot,args.amplicon_seq,cut_point, df_allele_around_cut,sgRNA,OUTPUT_DIRECTORY,MIN_FREQUENCY=args.min_frequency_alleles_around_cut_to_plot,MAX_N_ROWS=args.max_rows_alleles_around_cut_to_plot)

             print('==>Done.')

             if not args.keep_intermediate:
                 print('Removing Intermediate files...')

                 if args.fastq_r2!='':
                     files_to_remove=[processed_output_filename,flash_hist_filename,flash_histogram_filename,\
                                  flash_not_combined_1_filename,flash_not_combined_2_filename,\
                                  database_fasta_filename]
                 else:
                     files_to_remove=[processed_output_filename,database_fasta_filename]

                 if args.trim_sequences and args.fastq_r2!='':
                     files_to_remove+=[output_forward_paired_filename,output_reverse_paired_filename,\
                                                       output_forward_unpaired_filename,output_reverse_unpaired_filename]

                 if not args.dump:
                     files_to_remove+=[needle_output_filename]
                     if args.expected_hdr_amplicon_seq:
                         files_to_remove+=[needle_output_repair_filename]

                 if args.expected_hdr_amplicon_seq:
                     files_to_remove+=[database_repair_fasta_filename,]

                 if args.split_paired_end:
                     files_to_remove+=splitted_files_to_remove

                 if args.min_average_read_quality>0 or args.min_single_bp_quality>0:

                    if args.fastq_r2!='':
                             files_to_remove+=[args.fastq_r1,args.fastq_r2]
                    else:
                             files_to_remove+=[args.fastq_r1]

                 if sr_not_aligned.count():
                     files_to_remove+=[fasta_not_aligned_filename,database_rc_fasta_filename,needle_output_rc_filename]

                     if args.expected_hdr_amplicon_seq:
                            files_to_remove+=[database_repair_rc_fasta_filename,needle_output_repair_rc_filename]


                 for file_to_remove in files_to_remove:
                     try:
                             if os.path.islink(file_to_remove):
                                 os.unlink(file_to_remove)
                             else:
                                 os.remove(file_to_remove)
                     except:
                             warn('Skipping:%s' %file_to_remove)

             #write effect vectors as plain text files
             append = False
             print('Saving processed data...')
             def save_vector_to_file(vector,name):
                     np.savetxt(_jp('%s.txt' %name), np.vstack([(np.arange(len(vector))+1),vector]).T, fmt=['%d','%.18e'],delimiter='\t', newline='\n', header='amplicon position\teffect',footer='', comments='# ')

             # OUTPUT HEADER TO outputsummary.txt
             ref_donor_diffs=[i for i in xrange(len(args.expected_hdr_amplicon_seq)) if args.expected_hdr_amplicon_seq[i] != args.amplicon_seq[i]]
             filename=os.path.join(os.path.abspath(args.output_folder),args.name,'summary_of_editing_frequency.txt')
             fh_outfile = open(filename,"w")
             fh_outfile.write('Sample\tTotalReads\tMergedReads\tPercentMerged\tAlignedReads\tPercentAligned\tUnmodified\t%Unmodified\tCutsiteSubs\tNon-cutsiteSubs\t%CutsiteSubs\tNHEJ\t%NHEJ\t')

#-------------------------------------
             if os.path.isfile('outputsummary.txt'):
                 fh_summary = open('outputsummary.txt',"a")
                 append = True
             else:
                 fh_summary = open('outputsummary.txt',"w")
                 append = False
                 fh_summary.write('Sample\tTotalReads\tMergedReads\tPercentMerged\tAlignedReads\tPercentAligned\tUnmodified\t%Unmodified\tCutsiteSubs\tNon-cutsiteSubs\t%CutsiteSubs\tNHEJ\t%NHEJ\t')

             l_ref_donor = len(ref_donor_diffs);
             for i in range(len(ref_donor_diffs)):
                 if not args.main_site == i:
                     if not append:
                         fh_summary.write('Edit Site %d\tEdit %%\t' % (i+1) )
                     fh_outfile.write('Edit Site %d\tEdit %%\t' % (i+1) )
                 else: 
                     if not append:
                         fh_summary.write('Main HDR Site\tMain %\t' )
                     fh_outfile.write('Main HDR Site\tMain %\t' )
             if not append and (len(ref_donor_diffs) == 0):
                 fh_summary.write('\n')
             elif not append: 
                 fh_summary.write('All_HDR_Pos_Edited\tAll_Edit%\n')
             fh_outfile.write('All_HDR_Pos_Edited\tAll_Edit%\n')

             # SUMMARIZE
             old_denominator = float(N_INDELS + N_UNMODIFIED + N_CUTSUBS)
             percent_indel = 100*(N_INDELS / float(N_ALIGNED))
             percent_subs = 100*(N_CUTSUBS / float(N_ALIGNED))
             percent_unmod = 100*(N_UNMODIFIED / float(N_ALIGNED))
             # percent_aligned = 100 * (N_ALIGNED/float(N_READS_INPUT)) change to aligned over total merged
             percent_aligned = 100 * (N_ALIGNED / float(N_READS_AFTER_PREPROCESSING))
             percent_merged =  100 * (N_READS_AFTER_PREPROCESSING / float(N_READS_INPUT))
             
             fh_outfile.write('%s\t%d\t%d\t%.3f\t%d\t%.3f\t%d\t%.3f\t%d\t%d\t%.3f\t%d\t%.3f'  % (args.name,N_READS_INPUT,N_READS_AFTER_PREPROCESSING,percent_merged,N_ALIGNED,percent_aligned,N_UNMODIFIED,percent_unmod,N_CUTSUBS,seq_error,percent_subs,N_INDELS,percent_indel))
             fh_summary.write('%s\t%d\t%d\t%.3f\t%d\t%.3f\t%d\t%.3f\t%d\t%d\t%.3f\t%d\t%.3f'  % (args.name,N_READS_INPUT,N_READS_AFTER_PREPROCESSING,percent_merged,N_ALIGNED,percent_aligned,N_UNMODIFIED,percent_unmod,N_CUTSUBS,seq_error,percent_subs,N_INDELS,percent_indel))

             if args.expected_hdr_amplicon_seq:
                 last = len(ref_donor_diffs)
                 for i in range(len(ref_donor_diffs)):
                     percent_edit = 100 * (hdr_vector[i] / float(N_ALIGNED))
                     fh_outfile.write( ('\t%d\t%.3f' % (hdr_vector[i],percent_edit)))
                     fh_summary.write( ('\t%d\t%.3f' % (hdr_vector[i],percent_edit)))
                 # any site edited
                 percent_all_edit = 100 * (hdr_vector[last] / float(N_ALIGNED))
                 fh_outfile.write('\t%d\t%.3f' % (hdr_vector[last],percent_all_edit))
                 fh_summary.write('\t%d\t%.3f' % (hdr_vector[last],percent_all_edit))
             fh_outfile.write('\n')
             fh_summary.write('\n')
             fh_summary.close()
             fh_outfile.close()

             #write statistics
             if not args.skip_aln:
                 with open(_jp('Mapping_statistics.txt'),'w+') as outfile:
                     outfile.write('READS IN INPUTS:%d\nREADS AFTER PREPROCESSING:%d\nREADS ALIGNED:%d' % (N_READS_INPUT,N_READS_AFTER_PREPROCESSING,N_ALIGNED))

             with open(_jp('Quantification_of_editing_frequency.txt'),'w+') as outfile:
                     outfile.write(
                     ('Quantification of editing frequency\tUnmodified\tNHEJ\tHDR\n\t%d reads\t%d reads\t%d reads\n'  % (N_UNMODIFIED,N_INDELS,N_REPAIRED))\
                     +('\t- NHEJ:%d reads (%d reads with insertions, %d reads with deletions, %d reads with substitutions)\n' % (N_INDELS, np.sum(df_needle_alignment.ix[df_needle_alignment.NHEJ,'n_inserted']>0),np.sum(df_needle_alignment.ix[df_needle_alignment.NHEJ,'n_deleted']>0),np.sum(df_needle_alignment.ix[df_needle_alignment.NHEJ,'n_mutated']>0)))\
                     +('\t- HDR:%d reads (%d reads with insertions, %d reads with deletions, %d reads with substitutions)\n' % (N_REPAIRED, np.sum(df_needle_alignment.ix[df_needle_alignment.HDR,'n_inserted']>0),np.sum(df_needle_alignment.ix[df_needle_alignment.HDR,'n_deleted']>0),np.sum(df_needle_alignment.ix[df_needle_alignment.HDR,'n_mutated']>0)))\
                     +('Total Aligned:%d reads ' % N_ALIGNED))
             print('Done counting and outputting summary...')


#+('\t- Mixed HDR-NHEJ:%d reads (%d reads with insertions, %d reads with deletions, %d reads with substitutions)\n\n' % (N_MIXED_HDR_NHEJ, np.sum(df_needle_alignment.ix[df_needle_alignment.MIXED,'n_inserted']>0),np.sum(df_needle_alignment.ix[df_needle_alignment.MIXED,'n_deleted']>0),np.sum(df_needle_alignment.ix[df_needle_alignment.MIXED,'n_mutated']>0)))\


             #write alleles table
             df_alleles.ix[:,:'%Reads'].to_csv(_jp('Alleles_frequency_table.txt'),sep='\t',header=True,index=None)


             if PERFORM_FRAMESHIFT_ANALYSIS:
                 with open(_jp('Frameshift_analysis.txt'),'w+') as outfile:
                         outfile.write('Frameshift analysis:\n\tNoncoding mutation:%d reads\n\tIn-frame mutation:%d reads\n\tFrameshift mutation:%d reads\n' %(NON_INDELS_NON_FRAMESHIFT, MODIFIED_NON_FRAMESHIFT ,MODIFIED_FRAMESHIFT))

                 with open(_jp('Splice_sites_analysis.txt'),'w+') as outfile:
                         outfile.write('Splice sites analysis:\n\tUnmodified:%d reads\n\tPotential splice sites modified:%d reads\n' %(df_needle_alignment.shape[0]- SPLICING_SITES_MODIFIED, SPLICING_SITES_MODIFIED))


                 save_vector_to_file(effect_vector_insertion_noncoding,'effect_vector_insertion_noncoding')
                 save_vector_to_file(effect_vector_deletion_noncoding,'effect_vector_deletion_noncoding')
                 save_vector_to_file(effect_vector_mutation_noncoding,'effect_vector_substitution_noncoding')


             save_vector_to_file(effect_vector_insertion,'effect_vector_insertion_NHEJ')
             save_vector_to_file(effect_vector_deletion,'effect_vector_deletion_NHEJ')
             save_vector_to_file(effect_vector_mutation,'effect_vector_substitution_NHEJ')
             save_vector_to_file(effect_vector_combined,'effect_vector_combined')

             save_vector_to_file(avg_vector_ins_all,'position_dependent_vector_avg_insertion_size')
             save_vector_to_file(avg_vector_del_all,'position_dependent_vector_avg_deletion_size')


             pd.DataFrame(np.vstack([hlengths,hdensity]).T,columns=['indel_size','fq']).to_csv(_jp('indel_histogram.txt'),index=None,sep='\t')
             pd.DataFrame(np.vstack([x_bins_ins[:-1],y_values_ins]).T,columns=['ins_size','fq']).to_csv(_jp('insertion_histogram.txt'),index=None,sep='\t')
             pd.DataFrame(np.vstack([-x_bins_del[:-1],y_values_del]).T,columns=['del_size','fq']).to_csv(_jp('deletion_histogram.txt'),index=None,sep='\t')
             pd.DataFrame(np.vstack([x_bins_mut[:-1],y_values_mut]).T,columns=['sub_size','fq']).to_csv(_jp('substitution_histogram.txt'),index=None,sep='\t')

             if args.expected_hdr_amplicon_seq:
                 save_vector_to_file(effect_vector_insertion_mixed,'effect_vector_insertion_mixed_HDR_NHEJ')
                 save_vector_to_file(effect_vector_deletion_mixed,'effect_vector_deletion_mixed_HDR_NHEJ')
                 save_vector_to_file(effect_vector_mutation_mixed,'effect_vector_substitution_mixed_HDR_NHEJ')
                 save_vector_to_file(effect_vector_insertion_hdr,'effect_vector_insertion_HDR')
                 save_vector_to_file(effect_vector_deletion_hdr,'effect_vector_deletion_HDR')
                 save_vector_to_file(effect_vector_mutation_hdr,'effect_vector_substitution_HDR')


             if cut_points:
                cp.dump(sgRNA_intervals, open( _jp('sgRNA_intervals.pickle'), 'wb' ) )

             if sgRNA_intervals:
                cp.dump( cut_points, open( _jp('cut_points.pickle'), 'wb' ) )

             if offset_plots.any():
                  cp.dump(offset_plots,open( _jp('offset_plots.pickle'), 'wb' ) )

             if args.dump:
                 info('Dumping all the processed data...')
                 np.savez(_jp('effect_vector_insertion_NHEJ'),effect_vector_insertion)
                 np.savez(_jp('effect_vector_deletion_NHEJ'),effect_vector_deletion)
                 np.savez(_jp('effect_vector_substitution_NHEJ'),effect_vector_mutation)

                 np.savez(_jp('effect_vector_combined'),effect_vector_combined)

                 np.savez(_jp('position_dependent_vector_avg_insertion_size'),avg_vector_ins_all)
                 np.savez(_jp('position_dependent_vector_avg_deletion_size'),avg_vector_del_all)

                 df_needle_alignment.to_pickle(_jp('processed_reads_dataframe.pickle'))

                 if args.expected_hdr_amplicon_seq:
                     np.savez(_jp('effect_vector_insertion_mixed_HDR_NHEJ'),effect_vector_insertion_mixed)
                     np.savez(_jp('effect_vector_deletion_mixed_HDR_NHEJ'),effect_vector_deletion_mixed)
                     np.savez(_jp('effect_vector_substitution_mixed_HDR_NHEJ'),effect_vector_mutation_mixed)


                     np.savez(_jp('effect_vector_insertion_HDR'),effect_vector_insertion_hdr)
                     np.savez(_jp('effect_vector_deletion_HDR'),effect_vector_deletion_hdr)
                     np.savez(_jp('effect_vector_substitution_HDR'),effect_vector_mutation_hdr)

             print('==>All Done!')
             print('''
                      ))
                    ((
                      ))
                    _____
                 C\|~~~~~|
                    \___/
                ''')

             sys.exit(0)


    except NTException as e:
         output_error(e)
         error('NTException Alphabet error, please check your input.\n\nERROR: %s' % e)
         sys.exit(1)
    except SgRNASequenceException as e:
         output_error(e)
         error('SgRNASeqeuenceException sgRNA error, please check your input.\n\nERROR: %s' % e)
         sys.exit(2)
    except DonorSequenceException as e:
         output_error(e)
         error('Problem with the expected hdr amplicon sequence parameter, please check your input.\n\nERROR: %s' % e)
         sys.exit(3)
    except TrimmomaticException as e:
         output_error(e)
         error('Trimming error, please check your input.\n\nERROR: %s' % e)
         sys.exit(4)
    except ReorientException as e:
         output_error(e)
         error('Reorienting error.\n\nERROR: %s' % e)
         sys.exit(5)
    except FlashException as e:
         output_error(e)
         error('Merging error, please check your input.\n\nERROR: %s' % e)
         sys.exit(5)
    except NeedleException as e:
         output_error(e)
         error('Alignment error, please check your input.\n\nERROR: %s' % e)
         sys.exit(6)
    except NoReadsAlignedException as e:
         output_error(e)
         error('Alignment error, please check your input.\n\nERROR: %s' % e)
         sys.exit(7)
    except AmpliconEqualDonorException as e:
         output_error(e)
         error('Problem with the expected hdr amplicon sequence parameter, please check your input.\n\nERROR: %s' % e)
         sys.exit(8)
    except CoreDonorSequenceNotContainedException as e:
         output_error(e)
         error('Donor sequence error, please check your input.\n\nERROR: %s' % e)
         sys.exit(9)
    except CoreDonorSequenceNotUniqueException as e:
         output_error(e)
         error('Donor sequence error, please check your input.\n\nERROR: %s' % e)
         sys.exit(10)
    except ExonSequenceException as e:
         output_error(e)
         error('Coding sequence error, please check your input.\n\nERROR: %s' % e)
         sys.exit(11)
    except DuplicateSequenceIdException as e:
         output_error(e)
         error('Fastq file error, please check your input.\n\nERROR: %s' % e)
         sys.exit(12)
    except NoReadsAfterQualityFiltering as e:
         output_error(e)
         error('Filtering error, please check your input.\n\nERROR: %s' % e)
         sys.exit(13)
    except Exception as e:
         output_error(e)
         error('%s Unknown :-( error, please check your input.\n\nERROR: %s' % (args.name,e))
         sys.exit(-1)
