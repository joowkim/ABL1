import os
import subprocess
import numpy
import pandas

def onTargetReads(bam_file,amplicons_target_region_bed_file):
        cmd = '''samtools view -F 0x40 -c %s''' %(bam_file)
        mapped_reads = 'na'
        try:
            mapped_reads = int([x for x in subprocess.check_output(cmd, shell=True).split('\n') if x.strip() != ''][0])
        except:
            pass
        cmd = '''samtools view -L %s -c %s''' % (amplicons_target_region_bed_file,bam_file)
        mapped_to_target = int([x for x in subprocess.check_output(cmd, shell=True).split('\n') if x.strip() != ''][0])
        on_target = 'na'
        try:
            on_target = float(mapped_to_target) * 100 / mapped_reads
        except:
            pass
        return on_target,mapped_reads

def get_coverage(amplicons_target_region_bed_file,
                 bam_file,
                 minReadoverlap,
                 minAmpliconCov,
                 minSampleCov,
                 minReadUniform, 
                 output_coverage_file_tsv):

        cmd = '''bedtools coverage -a %s -b %s -counts -f %s > %s''' % (amplicons_target_region_bed_file,
                                                                        bam_file,
                                                                        str(minReadoverlap),
                                                                        output_coverage_file_tsv)
        subprocess.call(cmd, shell=True)
        Strg = '\\t'.join(
            ['Temp' + str(ix) if not ix in [3, 8] else ('target' if ix == 3 else 'coverage' + '\\n') for ix in
             range(9)])
        cmd = """sed -i '1s/^/%s/' %s""" % (Strg, output_coverage_file_tsv)
        subprocess.call(cmd, shell=True)
        dfcov = pandas.read_csv(output_coverage_file_tsv, sep='\t')

        mean_depth = numpy.mean(list(dfcov['coverage']))
        QC = 'Pass' if float(mean_depth) >= minSampleCov else 'Fail'

        dfcov['amplicon:coverage'] = dfcov['target'].astype(str)+':'+dfcov['coverage'].astype(str)
        dfcov['qc'] = dfcov['amplicon:coverage'].apply(lambda x: 'Pass' if int(x.split(':')[1]) >= minAmpliconCov else 'Fail')
        dfcov = dfcov.drop(columns=['amplicon:coverage'])
        df = dfcov.loc[dfcov['coverage'] >= minReadUniform*mean_depth]
        Uniformity = len(df) * 100.0 / len(dfcov)

        OnTarget,mapped_reads = onTargetReads(bam_file, amplicons_target_region_bed_file)

        Amplicons   = dfcov.drop(columns=[x for x in dfcov.columns if x.startswith('Temp')])
        Columns = ['mapped_reads', 'on_target', 'mean_depth', 'uniformity', 'qc']
        SummaryCov = [int(mapped_reads),round(OnTarget, 1),int(mean_depth), round(Uniformity,1),QC]
        SummaryCov = pandas.DataFrame([SummaryCov], columns=Columns)

if __name__ == "__main__":
  minAmpliconCov = 300
  minSampleCov = 100
  minReadUniform = 0.2
  minReadoverlap = 0.6
  amplicons_target_region_bed_file = ''
  bam_file = '' 
  output_coverage_file_tsv = ''

  get_coverage(amplicons_target_region_bed_file,
                 bam_file,
                 minReadoverlap,
                 minAmpliconCov,
                 minSampleCov,
                 minReadUniform,
                 output_coverage_file_tsv) 
