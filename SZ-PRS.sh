# SZ risk scores from the PGC
# 77 SNPs

plink --bfile /data/swe_gwas/ABZ/imputed_data/CNP/call_genos/cnp_callgenos_merged --score prs_77snp_input.txt sum --out SZ-PRS-77SNP-cnp
plink --bfile /data/swe_gwas/ABZ/imputed_data/OMNI/call_genos/omni_genos_merged --score prs_77snp_input_swe.txt sum --out SZ-PRS-77SNP-swe
