result_file=Tos_cis_eQTLs_1gPCs_quantile_10ePCs_lab_sex
awk -v bash_value="0.05" '{ if ($6 <= bash_value) print }' ${result_file} > file2.txt
