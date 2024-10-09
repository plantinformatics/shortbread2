outpdf="${genotypes/.vcf.gz/}-missing.pdf"
module load R
Rscript -e "\\n
     if(!require('vcfR'))\\n
         install.packages('vcfR')\\n

     if(!require('tidyverse'))\\n
         install.packages('tidyverse')\\n

     my_vcf <- read.vcfR('${genotypes}')\\n
     missingness <- extract.gt(my_vcf, elem = 'GT', which = 'all') \\%>\\%
               as.data.frame() \\%>\\%
               mutate(is_missing = ifelse(GT == './.', TRUE, FALSE))\\n

     missingness_per_sample <- missingness \\%>\\%
                           group_by(SAMPLE_ID) \\%>\\%
                           summarize(percent_missing = sum(is_missing) / n() * 100)\\n
     pdf('\${outpdf}',width=12,height=12)\\n
     ggplot(missingness_per_sample, aes(x = percent_missing)) +
     geom_histogram(binwidth = 1) +
     labs(x = 'Percentage of Missing Genotypes', y = 'Number of Samples') +
     theme_bw()\\n
     dev.off()\\n"