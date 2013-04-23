#!/usr/bin/env python
# File created on 23 Apr 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from cogent.maths.stats.test import pearson, spearman
from biom.parse import parse_biom_table
from qiime.util import parse_command_line_parameters, make_option

correlation_fs = {'spearman':spearman, 'pearson':pearson}

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = []
script_info['script_usage'].append(("","","%prog -i otu_table1.biom,otu_table2.biom -o pearson_correlated_otus.txt -m pearson"))
script_info['script_usage'].append(("","","%prog -i otu_table1.biom,otu_table2.biom -o spearman_correlated_otus.txt -m spearman"))
script_info['output_description']= ""

script_info['required_options'] = [
 make_option('-i','--input_biom_fps',type="existing_filepaths",help='the input filepaths'),
 make_option('-o','--output_fp',type="new_filepath",help='the output filepath'),
]

script_info['optional_options'] = [
 make_option('-m','--method',type="choice",choices=correlation_fs.keys(),
     help='the method to use for computing correlation [default: %default]', 
     default='spearman'),
]
script_info['version'] = __version__

def identify_correlated_observations(table1, table2, correlation_f=spearman):
    result = []
    shared_sample_ids = list(set(table1.SampleIds) & set(table2.SampleIds))
    
    def filter_sam_f(values,id,md):
        return id in shared_sample_ids
    def filter_obs_f(values,id,md):
        return values.sum() > 0
        
    table1_filt = table1.filterSamples(filter_sam_f).filterObservations(filter_obs_f).sortSampleOrder(shared_sample_ids)
    table2_filt = table2.filterSamples(filter_sam_f).filterObservations(filter_obs_f).sortSampleOrder(shared_sample_ids)
    
    for v1, i1, md1 in table1_filt.iterObservations():
        current_data = []
        for v2, i2, md2 in table2_filt.iterObservations():
            current_data.append(correlation_f(v1,v2))
        result.append(current_data)
        
    return result, table1_filt.ObservationIds, table2_filt.ObservationIds


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    data, row_ids, col_ids = identify_correlated_observations(
                              parse_biom_table(open(opts.input_biom_fps[0])), 
                              parse_biom_table(open(opts.input_biom_fps[1])),
                              correlation_f=correlation_fs[opts.method])
    output_f = open(opts.output_fp,'w')
    output_f.write('\t%s\n' % '\t'.join(col_ids))
    for row_id, datum in zip(row_ids, data):
        output_f.write('%s\t%s\n' % (row_id,'\t'.join(map(repr,datum))))
    output_f.close()
            

if __name__ == "__main__":
    main()