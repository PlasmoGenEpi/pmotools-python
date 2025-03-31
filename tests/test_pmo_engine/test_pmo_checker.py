#!/usr/bin/env python3

import os
import unittest
import json
import jsonschema
from pmotools.pmo_engine.pmo_checker import PMOChecker


class TestPMOChecker(unittest.TestCase):
    def setUp(self):
        self.working_dir = os.path.dirname(os.path.abspath(__file__))
        self.pmo_jsonschema_fnp = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "etc/portable_microhaplotype_object.schema.json")

        self.pmo_required_base_fields = ["experiment_info",
        "specimen_info",
        "sequencing_info",
        "panel_info",
        "target_info",
        "targeted_genomes",
        "microhaplotypes_info",
        "bioinformatics_methods_info",
        "bioinformatics_run_info",
        "microhaplotypes_detected",
        "pmo_header"]

        self.specimen_required_fields = [
                "specimen_name",
                "specimen_taxon_id",
                "host_taxon_id",
                "collection_date",
                "collection_country",
                "project_name"
            ]

    def test_pmo_checker_constructor(self):
        with open(self.pmo_jsonschema_fnp) as f: pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)

    def test_pmo_checker_check_for_required_base_fields(self):
        with open(self.pmo_jsonschema_fnp) as f: pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        pmo_test_object = {
            "experiment_info":[],
            "specimen_info":[],
            "sequencing_info":[],
            "panel_info":[],
            "target_info":[],
            "targeted_genomes":[],
            "microhaplotypes_info":{},
            "bioinformatics_methods_info":[],
            "bioinformatics_run_info":[],
            "microhaplotypes_detected":[],
            "pmo_header":{}
        }
        checker.check_for_required_base_fields(pmo_test_object)

    def test_pmo_checker_check_for_required_base_fields_fail(self):
        with open(self.pmo_jsonschema_fnp) as f: pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        pmo_test_object = {
            "experiment_info":[],
            "specimen_info":[],
            "sequencing_info":[],
            "panel_info":[],
            "target_info":[],
            "targeted_genome":[],
            "microhaplotypes_info":{},
            "bioinformatics_methods_infos":[],
            "bioinformatics_run_info":[],
            "microhaplotypes_detected":[],
            "pmo_headers":{}
        }
        self.assertRaises(Exception, checker.check_for_required_base_fields, pmo_test_object)

    def test_pmo_checker_get_required_fields_for_pmo_class(self):
        with open(self.pmo_jsonschema_fnp) as f: pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        self.assertEqual(checker.get_required_fields_for_pmo_class("SpecimenInfo"), self.specimen_required_fields)

    def test_pmo_checker_validate_pmo_json(self):
        with open(self.pmo_jsonschema_fnp) as f: pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        pmo_data = json.loads(
            '{"bioinformatics_methods_info":[{"additional_methods":[{"additional_argument":["--strictErrors","--illumina","--removeOneSampOnlyOneOffHaps","--excludeCommonlyLowFreqHaplotypes","--excludeLowFreqOneOffs","--rescueExcludedOneOffLowFreqHaplotypes"],"program":"SeekDeep processClusters","program_description":"Compare across samples for each target to create population level identifiers and do post artifact cleanup","program_version":"v2.6.5"}],"bioinformatics_method_name":"SeekDeep","demultiplexing_method":{"program":"SeekDeep extractorPairedEnd","program_description":"Takes raw paired-end reads and demultiplexes on primers and does QC filtering","program_version":"v2.6.5"},"denoising_method":{"additional_argument":["--illumina","--qualThres 25,20","--trimFront 1","--trimBack 1"],"program":"SeekDeep qluster","program_description":"Takes sequences per sample per target and clusters them","program_version":"v2.6.5"}}],"bioinformatics_run_info":[{"bioinformatics_methods_id":0,"bioinformatics_run_name":"Mozambique2018-SeekDeep","run_date":"2022-04-01"}],"experiment_info":[{"experiment_sample_name":"8025874217","panel_id":0,"sequencing_info_id":0,"sequencing_prep_plate_info":{"plate_col":12,"plate_name":"8","plate_row":"C"},"specimen_id":0}],"microhaplotypes_detected":[{"bioinformatics_run_id":0,"experiment_samples":[{"experiment_sample_id":0,"target_results":[{"haps":[{"mhap_id":0,"reads":34463},{"mhap_id":1,"reads":11600}],"mhaps_target_id":0}]}]}],"microhaplotypes_info":{"targets":[{"microhaplotypes":[{"seq":"AACTTTTTTTATTTTTTTTGTCAATAGATAAATGATCAATATTTTCTATATTTAATCTATCAAGTATTTTTATATATCTATTATTTCTTTCTTCGATGGATAAATTATATGAATCAATATCCTTTCTTTCATCAACAAACTTTTTTATTGTTAACTCCATTTTTTTATTTA"},{"seq":"AACTTTTTTTATTTTTTTTGTCAATAGATAAATGATCAATATTTTCTATATTTAATCTATCAAGTATTTTTATATATCTATTATTTCTTTCTTCGATGGATAAATTATAAGAATCAATATCCTTTCTTTCATCAACAAACTTTTTTATTGTTAACTCCATTTTTTTATTTA"}],"target_id":0}]},"panel_info":[{"panel_name":"heomev1","reactions":[{"panel_targets":[0],"reaction_name":"full"}]}],"pmo_header":{"creation_date":"2025-03-31","generation_method":{"program_name":"elucidator combingallintopmojson","program_version":"1.1.1"},"pmo_version":"v1.0.0"},"read_counts_by_stage":[{"bioinformatics_run_id":0,"read_counts_by_experimental_sample_by_stage":[{"experiment_sample_id":0,"read_counts_for_targets":[{"stages":[{"read_count":46063,"stage":"final"}],"target_id":0}],"total_raw_count":6858739}]}],"sequencing_info":[{"library_kit":"TruSeq i5/i7 barcode primers","library_layout":"paired-end","library_screen":"40 uL reaction containing 10 uL of bead purified digested product, 18uL of nuclease-free water, 8uL of 5X secondary PCR master mix, and 5 uL of 10 uM TruSeq i5/i7 barcode primers","library_selection":"RANDOM","library_source":"GENOMIC","library_strategy":"WGS","nucl_acid_amp":"https://www.paragongenomics.com/targeted-sequencing/amplicon-sequencing/cleanplex-ngs-amplicon-sequencing/","nucl_acid_amp_date":"2019-07-15","nucl_acid_ext":"https://www.paragongenomics.com/targeted-sequencing/amplicon-sequencing/cleanplex-ngs-amplicon-sequencing/","nucl_acid_ext_date":"2019-07-15","pcr_cond":"10 min at 95C, 13 cycles for high density samples (or 15 cycles for low density samples) of 15 sec at 98C and 75 sec at 60C","seq_center":"UCSF","seq_date":"2019-07-15","seq_instrument_model":"NextSeq 550 instrument","seq_platform":"ILLUMINA","sequencing_info_name":"Mozambique2018"}],"specimen_info":[{"collection_country":"Mozambique","collection_date":"2018-06-07","collector_chief_scientist":"Greenhouse, Bryan","geo_admin3":"Inhassoro","host_taxon_id":9606,"lat_lon":"-21.5535,35.1819","parasite_density_info":[{"density":477719.34375,"method":"qpcr"}],"plate_info":{"plate_col":12,"plate_name":"8","plate_row":"C"},"project_name":"MOZ2018","specimen_collect_device":"dried blood spot","specimen_name":"8025874217","specimen_store_loc":"UCSF Greenhouse Lab","specimen_taxon_id":[5833]}],"target_info":[{"forward_primers":[{"location":{"chrom":"Pf3D7_01_v3","end":145449,"genome_id":0,"start":145416,"strand":"+"},"seq":"TGTTCGATATGTTTAAATATATGATTCTCGAAA"}],"gene_name":"PF3D7_0103300","insert_location":{"chrom":"Pf3D7_01_v3","end":145622,"genome_id":0,"ref_seq":"AAACTTTTTTTATTTTTTTTGTCAATAGATAAATGATCAATATTTTCTATATTTAATCTATCAAGTATTTTTATATATCTATTATTTCTTTCTTCGATGGATAAATTATAAGAATCAATATCCTTTCTTTCATCAACAAACTTTTTTATTGTTAACTCCATTTTTTTATTTAA","start":145449,"strand":"+"},"reverse_primers":[{"location":{"chrom":"Pf3D7_01_v3","end":145449,"genome_id":0,"start":145416,"strand":"+"},"seq":"CCAATATGTCAAGGTATATTAAAGTATGGTATC"}],"target_attributes":["Included","Microhaplotype"],"target_name":"t1"}],"targeted_genomes":[{"chromosomes":["Pf3D7_01_v3","Pf3D7_02_v3","Pf3D7_03_v3","Pf3D7_04_v3","Pf3D7_05_v3","Pf3D7_06_v3","Pf3D7_07_v3","Pf3D7_08_v3","Pf3D7_09_v3","Pf3D7_10_v3","Pf3D7_11_v3","Pf3D7_12_v3","Pf3D7_13_v3","Pf3D7_14_v3","Pf3D7_API_v3","Pf3D7_MIT_v3"],"genome_version":"2020-09-01","gff_url":"https://plasmodb.org/common/downloads/release-65/Pfalciparum3D7/gff/data/PlasmoDB-65_Pfalciparum3D7.gff","name":"3D7","taxon_id":5833,"url":"https://plasmodb.org/common/downloads/release-65/Pfalciparum3D7/fasta/data/PlasmoDB-65_Pfalciparum3D7_Genome.fasta"}]}')
        checker.validate_pmo_json(pmo_data)

    def test_pmo_checker_validate_pmo_json_fail(self):
        with open(self.pmo_jsonschema_fnp) as f: pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        pmo_data = json.loads(
            '{"bioinformatics_methods_infos":[{"additional_methods":[{"additional_argument":["--strictErrors","--illumina","--removeOneSampOnlyOneOffHaps","--excludeCommonlyLowFreqHaplotypes","--excludeLowFreqOneOffs","--rescueExcludedOneOffLowFreqHaplotypes"],"program":"SeekDeep processClusters","program_description":"Compare across samples for each target to create population level identifiers and do post artifact cleanup","program_version":"v2.6.5"}],"bioinformatics_method_name":"SeekDeep","demultiplexing_method":{"program":"SeekDeep extractorPairedEnd","program_description":"Takes raw paired-end reads and demultiplexes on primers and does QC filtering","program_version":"v2.6.5"},"denoising_method":{"additional_argument":["--illumina","--qualThres 25,20","--trimFront 1","--trimBack 1"],"program":"SeekDeep qluster","program_description":"Takes sequences per sample per target and clusters them","program_version":"v2.6.5"}}],"bioinformatics_run_info":[{"bioinformatics_methods_id":0,"bioinformatics_run_name":"Mozambique2018-SeekDeep","run_date":"2022-04-01"}],"experiment_info":[{"experiment_sample_name":"8025874217","panel_id":0,"sequencing_info_id":0,"sequencing_prep_plate_info":{"plate_col":12,"plate_name":"8","plate_row":"C"},"specimen_id":0}],"microhaplotypes_detected":[{"bioinformatics_run_id":0,"experiment_samples":[{"experiment_sample_id":0,"target_results":[{"haps":[{"mhap_id":0,"reads":34463},{"mhap_id":1,"reads":11600}],"mhaps_target_id":0}]}]}],"microhaplotypes_info":{"targets":[{"microhaplotypes":[{"seq":"AACTTTTTTTATTTTTTTTGTCAATAGATAAATGATCAATATTTTCTATATTTAATCTATCAAGTATTTTTATATATCTATTATTTCTTTCTTCGATGGATAAATTATATGAATCAATATCCTTTCTTTCATCAACAAACTTTTTTATTGTTAACTCCATTTTTTTATTTA"},{"seq":"AACTTTTTTTATTTTTTTTGTCAATAGATAAATGATCAATATTTTCTATATTTAATCTATCAAGTATTTTTATATATCTATTATTTCTTTCTTCGATGGATAAATTATAAGAATCAATATCCTTTCTTTCATCAACAAACTTTTTTATTGTTAACTCCATTTTTTTATTTA"}],"target_id":0}]},"panel_info":[{"panel_name":"heomev1","reactions":[{"panel_targets":[0],"reaction_name":"full"}]}],"pmo_header":{"creation_date":"2025-03-31","generation_method":{"program_name":"elucidator combingallintopmojson","program_version":"1.1.1"},"pmo_version":"v1.0.0"},"read_counts_by_stage":[{"bioinformatics_run_id":0,"read_counts_by_experimental_sample_by_stage":[{"experiment_sample_id":0,"read_counts_for_targets":[{"stages":[{"read_count":46063,"stage":"final"}],"target_id":0}],"total_raw_count":6858739}]}],"sequencing_info":[{"library_kit":"TruSeq i5/i7 barcode primers","library_layout":"paired-end","library_screen":"40 uL reaction containing 10 uL of bead purified digested product, 18uL of nuclease-free water, 8uL of 5X secondary PCR master mix, and 5 uL of 10 uM TruSeq i5/i7 barcode primers","library_selection":"RANDOM","library_source":"GENOMIC","library_strategy":"WGS","nucl_acid_amp":"https://www.paragongenomics.com/targeted-sequencing/amplicon-sequencing/cleanplex-ngs-amplicon-sequencing/","nucl_acid_amp_date":"2019-07-15","nucl_acid_ext":"https://www.paragongenomics.com/targeted-sequencing/amplicon-sequencing/cleanplex-ngs-amplicon-sequencing/","nucl_acid_ext_date":"2019-07-15","pcr_cond":"10 min at 95C, 13 cycles for high density samples (or 15 cycles for low density samples) of 15 sec at 98C and 75 sec at 60C","seq_center":"UCSF","seq_date":"2019-07-15","seq_instrument_model":"NextSeq 550 instrument","seq_platform":"ILLUMINA","sequencing_info_name":"Mozambique2018"}],"specimen_info":[{"collection_country":"Mozambique","collection_date":"2018-06-07","collector_chief_scientist":"Greenhouse, Bryan","geo_admin3":"Inhassoro","host_taxon_id":9606,"lat_lon":"-21.5535,35.1819","parasite_density_info":[{"density":477719.34375,"method":"qpcr"}],"plate_info":{"plate_col":12,"plate_name":"8","plate_row":"C"},"project_name":"MOZ2018","specimen_collect_device":"dried blood spot","specimen_name":"8025874217","specimen_store_loc":"UCSF Greenhouse Lab","specimen_taxon_id":[5833]}],"target_info":[{"forward_primers":[{"location":{"chrom":"Pf3D7_01_v3","end":145449,"genome_id":0,"start":145416,"strand":"+"},"seq":"TGTTCGATATGTTTAAATATATGATTCTCGAAA"}],"gene_name":"PF3D7_0103300","insert_location":{"chrom":"Pf3D7_01_v3","end":145622,"genome_id":0,"ref_seq":"AAACTTTTTTTATTTTTTTTGTCAATAGATAAATGATCAATATTTTCTATATTTAATCTATCAAGTATTTTTATATATCTATTATTTCTTTCTTCGATGGATAAATTATAAGAATCAATATCCTTTCTTTCATCAACAAACTTTTTTATTGTTAACTCCATTTTTTTATTTAA","start":145449,"strand":"+"},"reverse_primers":[{"location":{"chrom":"Pf3D7_01_v3","end":145449,"genome_id":0,"start":145416,"strand":"+"},"seq":"CCAATATGTCAAGGTATATTAAAGTATGGTATC"}],"target_attributes":["Included","Microhaplotype"],"target_name":"t1"}],"targeted_genomes":[{"chromosomes":["Pf3D7_01_v3","Pf3D7_02_v3","Pf3D7_03_v3","Pf3D7_04_v3","Pf3D7_05_v3","Pf3D7_06_v3","Pf3D7_07_v3","Pf3D7_08_v3","Pf3D7_09_v3","Pf3D7_10_v3","Pf3D7_11_v3","Pf3D7_12_v3","Pf3D7_13_v3","Pf3D7_14_v3","Pf3D7_API_v3","Pf3D7_MIT_v3"],"genome_version":"2020-09-01","gff_url":"https://plasmodb.org/common/downloads/release-65/Pfalciparum3D7/gff/data/PlasmoDB-65_Pfalciparum3D7.gff","name":"3D7","taxon_id":5833,"url":"https://plasmodb.org/common/downloads/release-65/Pfalciparum3D7/fasta/data/PlasmoDB-65_Pfalciparum3D7_Genome.fasta"}]}')
        self.assertRaises(jsonschema.exceptions.ValidationError, checker.validate_pmo_json, pmo_data)


if __name__ == "__main__":
    unittest.main()
