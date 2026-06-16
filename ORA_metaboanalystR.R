library(MetaboAnalystR)
mSet <- InitDataObjects("conc", "pathora", FALSE, default.dpi = 300)

# 1. Read your data using standard R
temp_data <- read.csv('data/production/processed/plsda_metaboanalyst_ora.csv', header = TRUE)

# 2. Force row names to be unique (adds .1, .2, etc. to duplicates)
rownames(temp_data) <- make.unique(rownames(temp_data))

# 3. Save it to a temporary file
write.csv(temp_data, 'data/production/processed/temp_fixed_data.csv')

# 4. Try reading it with MetaboAnalystR again using the fixed file
mSet <- Read.TextData(mSet, 'data/production/processed/temp_fixed_data.csv', 'rowu', 'disc')


mSet <- CrossReferencing(mSet, 'Name')

mSet <- CreateMappingResultTable(mSet)
mSet <- SetKEGG.PathLib(mSet, 'hsa')
mSet <- CalculateOraScore(mSet, 'rbc', 'hyperg')
# -> 各 mechanism_class でサブセットして別々に実行することを推奨
