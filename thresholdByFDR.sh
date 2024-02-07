ls | grep -v '\.png$' | grep 'cis_eQTL' > all_result_files.txt

awk -v bash_value="0.05" '{ if ($6 <= bash_value) print }' ${result_file} > ${file}_0.05_FDR


while IFS= read -r line; do
    awk -v bash_value="0.05" '{ if ($6 <= bash_value) print }' $line > ${line}_0.05_FDR
done < all_result_files.txt
