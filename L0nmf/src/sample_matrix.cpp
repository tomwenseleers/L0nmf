#include <Rcpp.h>
#include <fstream>
#include <filesystem>
#include <vector>
#include <cstdint>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix slcn_read_block(std::string path, int n, IntegerVector channels,
                             IntegerVector rows) {
  const int nr=rows.size(),nc=channels.size();
  if(n<1||nr<1||nc<1)stop("Invalid matrix block dimensions.");
  for(int i=0;i<nr;++i)if(rows[i]<1||rows[i]>n)stop("Invalid matrix row.");
  for(int j=0;j<nc;++j)if(channels[j]<1)stop("Invalid matrix channel.");
  std::ifstream file(std::filesystem::u8path(path),std::ios::in|std::ios::binary);
  if(!file)stop("Cannot open sample matrix.");
  bool contiguous=true;
  for(int i=1;i<nr;++i)if(rows[i]!=rows[i-1]+1)contiguous=false;
  NumericMatrix out(nr,nc);
  std::vector<double> buffer(contiguous?0:n);
  const std::uint16_t endian=1;
  const bool little=*reinterpret_cast<const unsigned char*>(&endian)==1;
  for(int j=0;j<nc;++j) {
    const std::streamoff start=(static_cast<std::streamoff>(channels[j]-1)*n+
      (contiguous?rows[0]-1:0))*sizeof(double);
    const int count=contiguous?nr:n;
    double* dest=contiguous?out.begin()+static_cast<R_xlen_t>(j)*nr:buffer.data();
    file.seekg(start);file.read(reinterpret_cast<char*>(dest),static_cast<std::streamsize>(count)*sizeof(double));
    if(!file)stop("Truncated sample matrix.");
    if(!little)for(int i=0;i<count;++i) {
      unsigned char* p=reinterpret_cast<unsigned char*>(dest+i);
      for(int b=0;b<4;++b)std::swap(p[b],p[7-b]);
    }
    if(!contiguous)for(int i=0;i<nr;++i)out(i,j)=buffer[rows[i]-1];
  }
  return out;
}
