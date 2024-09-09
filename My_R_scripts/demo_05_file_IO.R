
# turn what to write out into dataframe
for_write_df = data.frame(each_seq_id, bin_subfolder_only_name, each_bin)

# write out into file
write.table(for_write_df, file_out_handle, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

