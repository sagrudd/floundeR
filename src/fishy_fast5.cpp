#include <Rcpp.h>
using namespace Rcpp;
#include<Rinternals.h>
#include <string.h>
using namespace std;
#include <zlib.h>
#include <stdio.h>
#include <cmath>

static std::string fastq_filename;
static int isZipped = 0;
static gzFile gzfile;
static FILE *file;
static int LONGEST_SEQ = 10000000;
static bool validFastq = true;
static size_t seqlength = 0;
static double seqqual = 0;

static vector<int> id_vect; 
static vector<int> length_vect; 
static vector<double> quality_vect; 
static std::vector<bool> qc_pass_vect;

struct Fastq_tag {
  char header[10000000];
  char sequence[10000000];
  char delim[10000000];
  char quality[10000000];
};

static Fastq_tag fq;



void init() {
  id_vect.clear();
  length_vect.clear();
  quality_vect.clear();
  qc_pass_vect.clear();
}



inline bool myfile_exists (const std::string& name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
}



inline int suffix_match(const char *query, const char *suffix) {
  return(strncmp(query + strlen(query) -
         strlen(suffix), suffix, strlen(suffix)));
}




inline int is_gzipped(std::string query)
{
  vector<string> list;
  list.push_back(".gzip");
  list.push_back(".gz");
  for( vector<string>::const_iterator it = list.begin(); it!=list.end(); ++it )
  {
    std::string suffix = *it;
    if (suffix_match(query.c_str(), suffix.c_str()) == 0)
    {
      return(1);
    }
  }
  return(0);
}


int has_next_fastq()
{
  if (isZipped == 1) {
    // Rcout << "hasNext(gzeof==" << gzeof(gzfile) << ")" << std::endl;
    if ((gzfile!=NULL) && (gzeof(gzfile)==0)) {
      return (1);
    }
  } else {
    if ((file!=NULL) && (feof(file)==0)) {
      return (1);
    }
  }
  return (0);
}

//' calculate mean Phred score from an ASCII encoded phred string
//' 
//' FASTQ and BAM store per base qualities as an ASCII string. Accessory methods
//' in e.g. ShortRead allow for a sum of the numeric encoded scores; this is not
//' corrected for the log/linear so scores are synthetically boosted
//' - this simple method performs mean on the character level data ...
//'
//' @param qstring is an ASCII encoded Phred quality score
//'
//' @return mean phred scaled q-value
//'
//' @examples
//' sequence_quality("ABCDEFGH")
//'
//' @export
// [[Rcpp::export]]
double sequence_quality(std::string qstring) 
{
  double running_total = 0.0;
  for (unsigned i=0; i<qstring.length(); ++i)
  {
    running_total += pow(
      (double) 10.0, (double) ((int) qstring.at(i) - 33) / -10.0);
  }
  double run_tot_mean =  -10 * log10( 
    running_total / qstring.length() );
  return(run_tot_mean);
}


int validate_fastq()
{
  if (has_next_fastq()==0) {
    Rcout << "Reading beyond end of file ..." << std::endl;
    return(-1); // run over the end of the file ...
  }
  
  // does header start with @ - if not there may be misalignment?
  char headdelim = fq.header[0];
  if (headdelim!='@')
  {
    Rcout << "Malformed fastq entry ~ Line1[@] delim not present" << std::endl;
    validFastq = false;
    return(0); // this is not a valid fastq
  }
  
  // is delim == "+"? - * is used to indicate runover end of file ...
  if (strcmp(fq.delim, "+")!=0)
  {
    Rcout << "Malformed fastq entry ~ Line3[+] not present" << std::endl;
    validFastq = false;
    return(0); // this is not a valid fastq
  }
  
  // is sequence length > 0
  seqlength = strlen(fq.sequence);
  size_t qlength = strlen(fq.quality);
  seqqual = sequence_quality(fq.quality);
  if (seqlength == 0) {
    Rcout << "Malformed fastq entry ~ length(sequence)==0" << std::endl;
    validFastq = false;
    return(0); // this is not a valid fastq
  }
  
  // does sequence length == quality length
  if (seqlength != qlength) {
    Rcout << "Malformed fastq entry ~ length(sequence)!=length(quality)" << std::endl;
    validFastq = false;
    return(0);
  }
  
  return (1);
}






int get_next_fastq()
{
  if (isZipped == 1)
  {
    gzgets(gzfile, fq.header, LONGEST_SEQ);
    gzgets(gzfile, fq.sequence, LONGEST_SEQ);
    gzgets(gzfile, fq.delim, LONGEST_SEQ);
    gzgets(gzfile, fq.quality, LONGEST_SEQ);
  }
  else
  {
    if (fgets(fq.header, LONGEST_SEQ, file) != NULL) {}
    if (fgets(fq.sequence, LONGEST_SEQ, file) != NULL) {}
    if (fgets(fq.delim, LONGEST_SEQ, file) != NULL) {}
    if (fgets(fq.quality, LONGEST_SEQ, file) != NULL) {}
  }
  // clip any trailing newlines
  fq.header[strcspn(fq.header, "\r\n")] = 0;
  fq.sequence[strcspn(fq.sequence, "\r\n")] = 0;
  fq.delim[strcspn(fq.delim, "\r\n")] = 0;
  fq.quality[strcspn(fq.quality, "\r\n")] = 0;
  return (validate_fastq());
}




//' Consume a FASTQ file and collate basic characteristic observations
//' 
//' This method is a fast C++ implementation of a FASTQ parser aimed to perform
//' a quick collation of summary_statistic information that may be used to make
//' decisions as to how FASTQ sequence collections are crafted.
//'
//' @return DataFrame of FASTQ based sequence observations
//'
//' @examples
//' fastq <- system.file("extdata", "example.fastq.gz", package = "floundeR")
//' fishy_fastq(fastq)
//'
//' @export
// [[Rcpp::export]]
DataFrame fishy_fastq(std::string fastq) {
  
  init();
  
  fastq_filename = fastq;
  
  // TEST (1) - DOES THE SPECIFIED FILE EXIST
  if (!myfile_exists(fastq_filename))
  {
    Rcout << "FastqFileNotFound" << std::endl;
    //return(NumericVector::get_na());
  }
  
  if (is_gzipped(fastq_filename)==1)
  {
    isZipped = 1;
    gzfile = gzopen(fastq_filename.c_str(), "r");
  } else {
    isZipped = 0;
    file = fopen(fastq_filename.c_str(), "r");
  }

  int counter = 0;
  
  while (get_next_fastq() == 1)
  {
    // Rcout << "(s/q)Lengths==" << seqlength << "/" << seqqual << std::endl;
    id_vect.push_back(counter); 
    length_vect.push_back(seqlength);
    quality_vect.push_back(seqqual);
    qc_pass_vect.push_back(true);
    counter += 1;
  }
  
  if (isZipped) {
    gzclose(gzfile);
  } else {
    fclose(file);
  }
  
  //NumericVector output = wrap(id_vect, length_vect);
  //output.attr("dim") = Dimension(id_vect.size(), 1);
  
  //NumericMatrix c = NumericMatrix(counter, 3);
  
  DataFrame df = DataFrame::create( 
    Named("id") = id_vect , 
    _["sequence_length_template"] = length_vect ,
    _["mean_qscore_template"] = quality_vect ,
    _["passes_filtering"] = qc_pass_vect);
  
  return df;
}
