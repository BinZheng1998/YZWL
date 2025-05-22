https://github.com/NBISweden/AGAT
（1）先将STRG基因从之前的文件中提取出来，然后用AGAT修复。 \
（2）合并（1）的修复后的STRG的gff文件以及从7b中liftoff过来的胡须gff文件，生成新的gff文件，并使用AGAT进行修复。 \
（3）20250522_agat_config.yaml 将check_identical_isoforms: false，check_all_level3_locations: false和create_l3_for_l2_orphan: false进行了修改。
