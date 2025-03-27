using ArgParse
using YRUtils
using JSON
using CSV
using DataFrames

TEMPLATES = Dict("atac_pipeline_wdl" => "/data/softwares/encode_pipeline/atac-seq-pipeline_v2.2.3/atac.wdl",
    "atac_input_json" => "/data/softwares/encode_pipeline/atac-seq-pipeline_v2.2.3/json_templates/atac_input.json",
    "backend_conf" => "/data/softwares/encode_pipeline/caper/local.conf")

s = ArgParseSettings(description="Run ENCODE ATAC-Seq pipeline.")
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
    help = "'atac' for ATAC-Seq or 'dnase' for DNase-Seq"
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
    help = "p-value threshold for peak calling"
    arg_type = Float64
    action = "store_arg"
    default = 0.01

    "--smooth_win"
    help = "Size of smoothing window for MACS2"
    arg_type = Int
    action = "store_arg"
    default = 150

    "--enable_idr"
    help = "Enable IDR"
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

    "--multimapping"
    help = "The number of multi-mapping reads allowed"
    arg_type = Int
    action = "store_arg"
    default = 4

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

# run pipeline in either paired-end or single-end mode
@info "start running ENCODE ATAC-Seq pipeline over all samples ..."
if files_read_type == "paired"
    input_json = open(TEMPLATES["atac_input_json"], "r") do io
        JSON.parse(io)
    end

    input_json["atac.pipeline_type"] = args["pipeline_type"][1]
    input_json["atac.align_only"] = args["align_only"]
    input_json["atac.true_rep_only"] = args["true_rep_only"]
    input_json["atac.genome_tsv"] = args["metadata_tsv"][1]
    input_json["atac.paired_end"] = true
    input_json["atac.multimapping"] = args["multimapping"]
    input_json["atac.mapq_thresh"] = args["mapq_thresh"]
    input_json["atac.no_dup_removal"] = args["no_dup_removal"]
    input_json["atac.pval_thresh"] = args["pval_thresh"]
    input_json["atac.smooth_win"] = args["smooth_win"]
    input_json["atac.enable_idr"] = args["enable_idr"]
    input_json["atac.idr_thresh"] = args["idr_thresh"]
    input_json["atac.filter_chrs"] = if (isempty(args["filter_chrs"]) || isempty(strip(args["filter_chrs"][1])))
        String[]
    else
        String.(strip.(split(args["filter_chrs"][1], ",")))
    end

    for (id, reps) in files_dict
        input_json["atac.title"] = id

        for (rep_num, rep) in reps
            input_json[string("atac.fastqs_", rep_num, "_R1")] = rep["R1"]
            input_json[string("atac.fastqs_", rep_num, "_R2")] = rep["R2"]
        end

        output_json_path = string(id, "_", basename(TEMPLATES["atac_input_json"]))
        open(output_json_path, "w") do io
            JSON.print(io, input_json)
        end

        cmd = `$caper_path run $(TEMPLATES["atac_pipeline_wdl"]) -c $(TEMPLATES["backend_conf"]) -i $output_json_path --max-concurrent-tasks $(args["num_jobs"]) --singularity --no-deepcopy`
        @info string("running ", cmd, " ...")
        run(cmd; wait=true)
        @info string("running ", cmd, " done!")
    end
elseif files_read_type == "single"
    input_json = open(TEMPLATES["atac_input_json"], "r") do io
        JSON.parse(io)
    end

    input_json["atac.pipeline_type"] = args["pipeline_type"][1]
    input_json["atac.align_only"] = args["align_only"]
    input_json["atac.true_rep_only"] = args["true_rep_only"]
    input_json["atac.genome_tsv"] = args["metadata_tsv"][1]
    input_json["atac.paired_end"] = false
    input_json["atac.multimapping"] = args["multimapping"]
    input_json["atac.mapq_thresh"] = args["mapq_thresh"]
    input_json["atac.no_dup_removal"] = args["no_dup_removal"]
    input_json["atac.pval_thresh"] = args["pval_thresh"]
    input_json["atac.smooth_win"] = args["smooth_win"]
    input_json["atac.enable_idr"] = args["enable_idr"]
    input_json["atac.idr_thresh"] = args["idr_thresh"]
    input_json["atac.filter_chrs"] = if (isempty(args["filter_chrs"]) || isempty(strip(args["filter_chrs"][1])))
        String[]
    else
        String.(strip.(split(args["filter_chrs"][1], ",")))
    end

    for (id, reps) in files_dict
        input_json["atac.title"] = id

        for (rep_num, rep) in reps
            input_json[string("atac.fastqs_", rep_num, "_R1")] = rep
        end

        output_json_path = string(id, "_", basename(TEMPLATES["atac_input_json"]))
        open(output_json_path, "w") do io
            JSON.print(io, input_json)
        end

        cmd = `$caper_path run $(TEMPLATES["atac_pipeline_wdl"]) -c $(TEMPLATES["backend_conf"]) -i $output_json_path --max-concurrent-tasks $(args["num_jobs"]) --singularity --no-deepcopy`
        @info string("running ", cmd, " ...")
        run(cmd; wait=true)
        @info string("running ", cmd, " done!")
    end
else
    @error "invalid read type (read type can be either 'paired' or 'single')"
end
@info "running ENCODE ATAC-Seq pipeline over all samples done!"

# collect results
@info "start collecting outputs ..."
croo_path = YRUtils.ShellUtils.find_cmd("croo"; return_nothing=false)
YRUtils.ShellUtils.cmd_valid(Cmd(string.([croo_path, "--version"])); return_false=false)

metadata_json_pathes = YRUtils.BaseUtils.list_files(pwd(), r"metadata\.json$"; recursive=true, full_name=true)
atac_metadata_json_pathes = metadata_json_pathes[.!isnothing.(match.(r"(?<!\.caper_tmp)/atac/", metadata_json_pathes))]
for atac_metadata_json_path in atac_metadata_json_pathes
    atac_metadata_json = open(atac_metadata_json_path, "r") do io
        JSON.parse(io)
    end

    title = atac_metadata_json["inputs"]["title"]
    if !isnothing(title) && !isempty(title)
        output_dir = joinpath(args["output_dir"][1], title)
    else
        @error "cannot find the title field in metadata.json"
    end

    mkpath(output_dir)

    croo_cmd = `$croo_path $atac_metadata_json_path --out-dir $output_dir --method copy`
    @info string("running ", croo_cmd, " ...")
    run(croo_cmd; wait=true)
    @info string("running ", croo_cmd, " done!")
    cp(atac_metadata_json_path, joinpath(output_dir, basename(atac_metadata_json_path)))
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
