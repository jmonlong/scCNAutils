library(dplyr)

df = read.table('cells-community-normal-natcom2020.tsv.gz', as.is=TRUE, header=TRUE, sep='\t')

cat(unique(df$sample), sep='\n')

conv = c('OPK322'='BT322_GSC',
         'OPK324S'='BT324_GSC',
         'OPK326S'='BT326_GSC',
         'OPK333B'='BT333_W',
         'OPK333S'='BT333_GSC',
         'OPK338B'='BT338_W',
         'OPK346B'='BT346_W',
         'OPK363'='BT363_W',
         'OPK363s'='BT363_GSC',
         'OPK364B'='BT364_W',
         'OPK368B'='BT368_W',
         'OPK368S'='BT368_GSC',
         'OPK389B'='BT389_W',
         'OPK390'='BT390_W',
         'OPK397'='BT397_W',
         'OPK400'='BT400_W',
         'OPK402B'='BT402_W',
         'OPK407'='BT407_W',
         'OPK409'='BT409_W')

df2 = df %>% mutate(sample=conv[sample], community.patient=gsub('OPK', 'BT', community.patient))

write.table(df2, file='cells-community-normal-natcom2020.tsv', row.names=FALSE, quote=FALSE, sep='\t')
