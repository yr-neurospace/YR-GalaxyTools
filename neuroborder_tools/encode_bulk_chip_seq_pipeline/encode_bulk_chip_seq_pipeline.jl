using ArgParse
using YRUtils
using JSON
using CSV
using DataFrames

TEMPLATES = Dict("chip_pipeline_wdl" => "/data/softwares/encode_pipeline/chip-seq-pipeline2_v2.2.2/chip.wdl",
    "chip_input_json" => "/data/softwares/encode_pipeline/chip-seq-pipeline2_v2.2.2/json_templates/chip_input.json",
    "backend_conf" => "/data/softwares/encode_pipeline/caper/local.conf")

s = ArgParseSettings(description="Run ENCODE ChIP-Seq pipeline.")
@add_arg_table s begin
    "--names"
    help = "input file names separated by white space"
    nargs = '+'
    arg_type = String
    action = "append_arg"
    required = true

    "--files"
    help = "input file paths separated by white space"
    nargs = '+'
    arg_type = String
    action = "append_arg"
    required = true

    "--filter_chrs"
    help = "Chromosome names, separated with comma, to be filtered out from a final BAM file (e.g. 'chrM,MT')"
    nargs = 1
    arg_type = String
    action = "store_arg"
    required = false

    "--treatment_ctl_id_pairs"
    help = "Specify treatment-control ID pairs (paired treatment and control IDs are separated with colon, and different treatment-control pairs are separated with comma, e.g. 'Test1:Input1,Test2:Input2')"
    nargs = 1
    arg_type = String
    action = "store_arg"
    required = true

    "--input_dir"
    help = "input directory (i.e. all FASTQ files must be within this directory)"
    nargs = 1
    arg_type = String
    action = "store_arg"
    required = true

    "--output_dir"
    help = "output directory"
    nargs = 1
    arg_type = String
    action = "store_arg"
    required = true

    "--metadata_tsv"
    help = "path to genome metadata TSV file"
    nargs = 1
    arg_type = String
    action = "store_arg"
    required = true

    "--pipeline_type"
    help = "'tf' for TF ChIP-Seq, 'histone' for hitone ChIP-Seq, or 'control' for only mapping control FASTQs only"
    nargs = 1
    arg_type = String
    action = "store_arg"
    required = true

    "--peak_caller"
    help = "'spp' (SPP) for 'tf' (TF ChIP-Seq) type, and 'macs2' (MACS2) for 'histone' (histone ChIP-Seq) type"
    nargs = 1
    arg_type = String
    action = "store_arg"
    required = true

    "--align_only"
    help = "Only align reads to the reference genome (i.e. peak calling and its downstream analyses will be disabled)"
    action = "store_true"

    "--true_rep_only"
    help = "Disable pseudo-replicate generation and all related analyses"
    action = "store_true"

    "--no_dup_removal"
    help = "Do not remove read duplicates"
    action = "store_true"

    "--pval_thresh"
    help = "p-value threshold for peak caller MACS2"
    arg_type = Float64
    action = "store_arg"
    default = 0.01

    "--fdr_thresh"
    help = "FDR threshold for peak-caller SPP"
    arg_type = Float64
    action = "store_arg"
    default = 0.01

    "--cap_num_peak"
    help = "Cap number of peaks called"
    arg_type = Int
    action = "store_arg"
    default = 500000

    "--always_use_pooled_ctl"
    help = "Always use a pooled control to compare with each replicate (or choose an appropriate control for each replicate)"
    action = "store_true"

    "--use_bowtie2_local_mode"
    help = "Use Bowtie2's local mode (the default is end-to-end mode)"
    action = "store_true"

    "--idr_thresh"
    help = "Threshold for IDR"
    arg_type = Float64
    action = "store_arg"
    default = 0.05

    "--mapq_thresh"
    help = "Threshold for mapped reads quality (samtools view -q)"
    arg_type = Int
    action = "store_arg"
    default = 30

    "--num_jobs"
    help = "the number of concurrent jobs"
    arg_type = Int
    action = "store_arg"
    default = 1
end

args = parse_args(s)
for k in keys(args)
    if typeof(args[k]) <: AbstractArray
        args[k] = YRUtils.BaseUtils.flatten_array(args[k])
    end
end
if args["pipeline_type"][1] ∉ ("tf", "histone", "control")
    @error "invalid pipeline type"
end
if args["peak_caller"][1] ∉ ("spp", "macs2")
    @error "invalid peak caller"
end
args["froms2tos"] = collect(zip(args["files"], joinpath.(args["input_dir"][1], args["names"])))

YRUtils.BaseUtils.file_transfer(first.(args["froms2tos"]), last.(args["froms2tos"]); ops=["soft", "hard", "copy"])

dict = YRUtils.BioUtils.auto_detect_fastq_read_type(last.(args["froms2tos"]); allow_mix=false)
files_dict = if dict["paired"]["status"] == "yes"
    dict["paired"]["dict"]
elseif dict["single"]["status"] == "yes"
    dict["single"]["dict"]
else
    @error "did not detect any paired-end or single-end files"
end
files_read_type = if dict["paired"]["status"] == "yes"
    "paired"
elseif dict["single"]["status"] == "yes"
    "single"
else
    @error "did not detect any paired-end or single-end files"
end

caper_path = YRUtils.ShellUtils.find_cmd("caper"; return_nothing=false)
YRUtils.ShellUtils.cmd_valid(Cmd(string.([caper_path, "--version"])); return_false=false)

metadata_df = CSV.read(args["metadata_tsv"][1], DataFrame; header=["key", "value"], delim="\t")
metadata_dict = Dict(zip(metadata_df.key, metadata_df.value))

if args["pipeline_type"][1] ∈ ("tf", "histone")
    treatment_ctl_id_pairs = split.(split(args["treatment_ctl_id_pairs"][1], ","), ":")
    if any(length.(treatment_ctl_id_pairs) .!= 2)
        @error "Treatment-Control ID pairs are invalid"
    end
    treatment_ctl_id_pairs = [[strip(treatment), strip(ctl)] for (treatment, ctl) in treatment_ctl_id_pairs]
    flattened_treatment_ctl_id_pairs = YRUtils.BaseUtils.flatten_array(treatment_ctl_id_pairs)
    if any(isempty.(flattened_treatment_ctl_id_pairs)) || (sort(flattened_treatment_ctl_id_pairs) != sort(collect(keys(files_dict))))
        @error "Treatment-Control ID pairs are invalid"
    end
    treatment_ctl_id_pairs = convert(Vector{Vector{String}}, treatment_ctl_id_pairs)
end

# run pipeline in either paired-end or single-end mode
@info "start running ENCODE ChIP-Seq pipeline over all samples ..."
if files_read_type == "paired"
    input_json = open(TEMPLATES["chip_input_json"], "r") do io
        JSON.parse(io)
    end

    input_json["chip.pipeline_type"] = args["pipeline_type"][1]
    input_json["chip.align_only"] = args["align_only"]
    input_json["chip.true_rep_only"] = args["true_rep_only"]
    input_json["chip.genome_tsv"] = args["metadata_tsv"][1]
    input_json["chip.paired_end"] = true
    input_json["chip.ctl_paired_end"] = true
    input_json["chip.mapq_thresh"] = args["mapq_thresh"]
    input_json["chip.no_dup_removal"] = args["no_dup_removal"]
    input_json["chip.pval_thresh"] = args["pval_thresh"]
    input_json["chip.idr_thresh"] = args["idr_thresh"]
    input_json["chip.peak_caller"] = args["peak_caller"][1]
    input_json["chip.use_bowtie2_local_mode"] = args["use_bowtie2_local_mode"]
    input_json["chip.always_use_pooled_ctl"] = args["always_use_pooled_ctl"]
    input_json["chip.cap_num_peak"] = args["cap_num_peak"]
    input_json["chip.fdr_thresh"] = args["fdr_thresh"]
    input_json["chip.filter_chrs"] = if (isempty(args["filter_chrs"]) || isempty(strip(args["filter_chrs"][1])))
        String[]
    else
        String.(strip.(split(args["filter_chrs"][1], ",")))
    end

    if args["pipeline_type"][1] ∈ ("tf", "histone")
        for (treatment, ctl) in treatment_ctl_id_pairs
            input_json["chip.title"] = treatment

            for (rep_num, rep) in files_dict[treatment]
                input_json[string("chip.fastqs_", rep_num, "_R1")] = rep["R1"]
                input_json[string("chip.fastqs_", rep_num, "_R2")] = rep["R2"]
            end

            for (rep_num, rep) in files_dict[ctl]
                input_json[string("chip.ctl_fastqs_", rep_num, "_R1")] = rep["R1"]
                input_json[string("chip.ctl_fastqs_", rep_num, "_R2")] = rep["R2"]
            end

            output_json_path = string(treatment, "_", basename(TEMPLATES["chip_input_json"]))
            open(output_json_path, "w") do io
                JSON.print(io, input_json)
            end

            cmd = `$caper_path run $(TEMPLATES["chip_pipeline_wdl"]) -c $(TEMPLATES["backend_conf"]) -i $output_json_path --max-concurrent-tasks $(args["num_jobs"]) --singularity --no-deepcopy`
            @info string("running ", cmd, " ...")
            run(cmd; wait=true)
            @info string("running ", cmd, " done!")
        end
    elseif args["pipeline_type"][1] ∈ ("control",)
        for (id, reps) in files_dict
            input_json["chip.title"] = id

            for (rep_num, rep) in reps
                input_json[string("chip.fastqs_", rep_num, "_R1")] = rep["R1"]
                input_json[string("chip.fastqs_", rep_num, "_R2")] = rep["R2"]
            end

            output_json_path = string(id, "_", basename(TEMPLATES["chip_input_json"]))
            open(output_json_path, "w") do io
                JSON.print(io, input_json)
            end

            cmd = `$caper_path run $(TEMPLATES["chip_pipeline_wdl"]) -c $(TEMPLATES["backend_conf"]) -i $output_json_path --max-concurrent-tasks $(args["num_jobs"]) --singularity --no-deepcopy`
            @info string("running ", cmd, " ...")
            run(cmd; wait=true)
            @info string("running ", cmd, " done!")
        end
    else
        @error "invalid pipeline type"
    end
elseif files_read_type == "single"
    input_json = open(TEMPLATES["chip_input_json"], "r") do io
        JSON.parse(io)
    end

    input_json["chip.pipeline_type"] = args["pipeline_type"][1]
    input_json["chip.align_only"] = args["align_only"]
    input_json["chip.true_rep_only"] = args["true_rep_only"]
    input_json["chip.genome_tsv"] = args["metadata_tsv"][1]
    input_json["chip.paired_end"] = false
    input_json["chip.ctl_paired_end"] = false
    input_json["chip.mapq_thresh"] = args["mapq_thresh"]
    input_json["chip.no_dup_removal"] = args["no_dup_removal"]
    input_json["chip.pval_thresh"] = args["pval_thresh"]
    input_json["chip.idr_thresh"] = args["idr_thresh"]
    input_json["chip.peak_caller"] = args["peak_caller"][1]
    input_json["chip.use_bowtie2_local_mode"] = args["use_bowtie2_local_mode"]
    input_json["chip.always_use_pooled_ctl"] = args["always_use_pooled_ctl"]
    input_json["chip.cap_num_peak"] = args["cap_num_peak"]
    input_json["chip.fdr_thresh"] = args["fdr_thresh"]
    input_json["chip.filter_chrs"] = if (isempty(args["filter_chrs"]) || isempty(strip(args["filter_chrs"][1])))
        String[]
    else
        String.(strip.(split(args["filter_chrs"][1], ",")))
    end

    if args["pipeline_type"][1] ∈ ("tf", "histone")
        for (treatment, ctl) in treatment_ctl_id_pairs
            input_json["chip.title"] = treatment

            for (rep_num, rep) in files_dict[treatment]
                input_json[string("chip.fastqs_", rep_num, "_R1")] = rep
            end

            for (rep_num, rep) in files_dict[ctl]
                input_json[string("chip.ctl_fastqs_", rep_num, "_R1")] = rep
            end

            output_json_path = string(treatment, "_", basename(TEMPLATES["chip_input_json"]))
            open(output_json_path, "w") do io
                JSON.print(io, input_json)
            end

            cmd = `$caper_path run $(TEMPLATES["chip_pipeline_wdl"]) -c $(TEMPLATES["backend_conf"]) -i $output_json_path --max-concurrent-tasks $(args["num_jobs"]) --singularity --no-deepcopy`
            @info string("running ", cmd, " ...")
            run(cmd; wait=true)
            @info string("running ", cmd, " done!")
        end
    elseif args["pipeline_type"][1] ∈ ("control",)
        for (id, reps) in files_dict
            input_json["chip.title"] = id

            for (rep_num, rep) in reps
                input_json[string("chip.fastqs_", rep_num, "_R1")] = rep
            end

            output_json_path = string(id, "_", basename(TEMPLATES["chip_input_json"]))
            open(output_json_path, "w") do io
                JSON.print(io, input_json)
            end

            cmd = `$caper_path run $(TEMPLATES["chip_pipeline_wdl"]) -c $(TEMPLATES["backend_conf"]) -i $output_json_path --max-concurrent-tasks $(args["num_jobs"]) --singularity --no-deepcopy`
            @info string("running ", cmd, " ...")
            run(cmd; wait=true)
            @info string("running ", cmd, " done!")
        end
    else
        @error "invalid pipeline type"
    end
else
    @error "invalid read type (read type can be either 'paired' or 'single')"
end
@info "running ENCODE ChIP-Seq pipeline over all samples done!"

# collect results
@info "start collecting outputs ..."
croo_path = YRUtils.ShellUtils.find_cmd("croo"; return_nothing=false)
YRUtils.ShellUtils.cmd_valid(Cmd(string.([croo_path, "--version"])); return_false=false)

metadata_json_pathes = YRUtils.BaseUtils.list_files(pwd(), r"metadata\.json$"; recursive=true, full_name=true)
chip_metadata_json_pathes = metadata_json_pathes[.!isnothing.(match.(r"(?<!\.caper_tmp)/chip/", metadata_json_pathes))]
for chip_metadata_json_path in chip_metadata_json_pathes
    chip_metadata_json = open(chip_metadata_json_path, "r") do io
        JSON.parse(io)
    end

    title = chip_metadata_json["inputs"]["title"]
    if !isnothing(title) && !isempty(title)
        output_dir = joinpath(args["output_dir"][1], title)
    else
        @error "cannot find the title field in metadata.json"
    end

    mkpath(output_dir)

    croo_cmd = `$croo_path $chip_metadata_json_path --out-dir $output_dir --method copy`
    @info string("running ", croo_cmd, " ...")
    run(croo_cmd; wait=true)
    @info string("running ", croo_cmd, " done!")
    cp(chip_metadata_json_path, joinpath(output_dir, basename(chip_metadata_json_path)))
end
@info "collecting outputs done!"

# packaging and compressing directories and files
@info "start packaging and compressing directories and files ..."
tar_path = YRUtils.ShellUtils.find_cmd("tar"; return_nothing=false)
YRUtils.ShellUtils.cmd_valid(Cmd(string.([tar_path, "--version"])); return_false=false)
pigz_path = YRUtils.ShellUtils.find_cmd("pigz"; return_nothing=false)
YRUtils.ShellUtils.cmd_valid(Cmd(string.([pigz_path, "--version"])); return_false=false)

cd(args["output_dir"][1])
for item in readdir(pwd(); join=false)
    if isdir(item)
        tar_gz_cmd = `$tar_path --use-compress-program=$pigz_path -cvf $(item).tgz $item`
        @info string("running ", tar_gz_cmd, " ...")
        run(tar_gz_cmd; wait=true)
        @info string("running ", tar_gz_cmd, " done!")
    end
end
@info "packaging and compressing directories and files done!"
