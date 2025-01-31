using ArgParse
using YRUtils
using JSON
using CSV
using DataFrames

TEMPLATES = Dict("rna_pipeline_wdl" => "/data/softwares/encode_pipeline/rna-seq-pipeline_v1.2.4/rna-seq-pipeline.wdl",
    "rna_pe_json_without_kallisto" => "/data/softwares/encode_pipeline/rna-seq-pipeline_v1.2.4/json_templates/rna_input_no_kallisto_pe.json",
    "rna_se_json_without_kallisto" => "/data/softwares/encode_pipeline/rna-seq-pipeline_v1.2.4/json_templates/rna_input_no_kallisto_se.json",
    "backend_conf" => "/data/softwares/encode_pipeline/caper/local.conf")

s = ArgParseSettings(description="Run ENCODE RNA-seq pipeline.")
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
    help = "path to genome metadata TSV file (in which STAR index, RSEM index, chromosome sizes, and transcript IDs to gene types table are given)"
    nargs = 1
    arg_type = String
    action = "store_arg"
    required = true

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
@info "start running ENCODE RNA-seq pipeline over all samples ..."
if files_read_type == "paired"
    input_json = open(TEMPLATES["rna_pe_json_without_kallisto"], "r") do io
        JSON.parse(io)
    end

    input_json["rna.align_index"] = metadata_dict["star_index"]
    input_json["rna.rsem_index"] = metadata_dict["rsem_index"]
    input_json["rna.chrom_sizes"] = metadata_dict["chrom_sizes"]
    input_json["rna.rna_qc_tr_id_to_gene_type_tsv"] = metadata_dict["tx2gene_type"]
    input_json["rna.endedness"] = files_read_type

    for (id, reps) in files_dict
        input_json["rna.bamroot"] = id

        fq_r1_files = Vector{String}[]
        fq_r2_files = Vector{String}[]
        for (_, reads) in reps
            push!(fq_r1_files, reads["R1"])
            push!(fq_r2_files, reads["R2"])
        end
        input_json["rna.fastqs_R1"] = fq_r1_files
        input_json["rna.fastqs_R2"] = fq_r2_files

        output_json_path = string(id, "_", basename(TEMPLATES["rna_pe_json_without_kallisto"]))
        open(output_json_path, "w") do io
            JSON.print(io, input_json)
        end

        cmd = `$caper_path run $(TEMPLATES["rna_pipeline_wdl"]) -c $(TEMPLATES["backend_conf"]) -i $output_json_path --max-concurrent-tasks $(args["num_jobs"]) --singularity --no-deepcopy`
        @info string("running ", cmd, " ...")
        run(cmd; wait=true)
        @info string("running ", cmd, " done!")
    end
elseif files_read_type == "single"
    input_json = open(TEMPLATES["rna_se_json_without_kallisto"], "r") do io
        JSON.parse(io)
    end

    input_json["rna.align_index"] = metadata_dict["star_index"]
    input_json["rna.rsem_index"] = metadata_dict["rsem_index"]
    input_json["rna.chrom_sizes"] = metadata_dict["chrom_sizes"]
    input_json["rna.rna_qc_tr_id_to_gene_type_tsv"] = metadata_dict["tx2gene_type"]
    input_json["rna.endedness"] = files_read_type

    for (id, reps) in files_dict
        input_json["rna.bamroot"] = id
        input_json["rna.fastqs_R1"] = [parts for (_, parts) in reps]

        output_json_path = string(id, "_", basename(TEMPLATES["rna_se_json_without_kallisto"]))
        open(output_json_path, "w") do io
            JSON.print(io, input_json)
        end

        caper_cmd = `$caper_path run $(TEMPLATES["rna_pipeline_wdl"]) -c $(TEMPLATES["backend_conf"]) -i $output_json_path --max-concurrent-tasks $(args["num_jobs"]) --singularity --no-deepcopy`
        @info string("running ", caper_cmd, " ...")
        run(caper_cmd; wait=true)
        @info string("running ", caper_cmd, " done!")
    end
else
    @error "invalid read type (read type can be either 'paired' or 'single')"
end
@info "running ENCODE RNA-seq pipeline over all samples done!"

# collect results
@info "start collecting outputs ..."
croo_path = YRUtils.ShellUtils.find_cmd("croo"; return_nothing=false)
YRUtils.ShellUtils.cmd_valid(Cmd(string.([croo_path, "--version"])); return_false=false)

metadata_json_pathes = YRUtils.BaseUtils.list_files(pwd(), r"metadata\.json$"; recursive=true, full_name=true)
rna_metadata_json_pathes = metadata_json_pathes[.!isnothing.(match.(r"/rna/", metadata_json_pathes))]
for rna_metadata_json_path in rna_metadata_json_pathes
    rna_metadata_json = open(rna_metadata_json_path, "r") do io
        JSON.parse(io)
    end

    # delete mad_qc-related fileds
    # croo cannot work with it correctly for some unknown reasons, if mad_qc was run
    # the only way that I can fix this bug is to delete the mad_qc-related fileds
    if haskey(rna_metadata_json, "calls") && haskey(rna_metadata_json["calls"], "rna.mad_qc")
        delete!(rna_metadata_json["calls"], "rna.mad_qc")
    end
    if haskey(rna_metadata_json, "outputs") && haskey(rna_metadata_json["outputs"], "rna.mad_qc.madQCmetrics")
        delete!(rna_metadata_json["outputs"], "rna.mad_qc.madQCmetrics")
    end
    if haskey(rna_metadata_json, "outputs") && haskey(rna_metadata_json["outputs"], "rna.mad_qc.python_log")
        delete!(rna_metadata_json["outputs"], "rna.mad_qc.python_log")
    end
    if haskey(rna_metadata_json, "outputs") && haskey(rna_metadata_json["outputs"], "rna.mad_qc.madQCplot")
        delete!(rna_metadata_json["outputs"], "rna.mad_qc.madQCplot")
    end
    if haskey(rna_metadata_json, "inputs") && haskey(rna_metadata_json["inputs"], "mad_qc_disk")
        delete!(rna_metadata_json["inputs"], "mad_qc_disk")
    end

    open(rna_metadata_json_path, "w") do io
        JSON.print(io, rna_metadata_json)
    end

    bamroot = rna_metadata_json["inputs"]["bamroot"]
    if !isnothing(bamroot) && !isempty(bamroot)
        output_dir = joinpath(args["output_dir"][1], bamroot)
    else
        @error "cannot find the bamroot field in metadata.json"
    end

    mkpath(output_dir)

    croo_cmd = `$croo_path $rna_metadata_json_path --out-dir $output_dir --method copy`
    @info string("running ", croo_cmd, " ...")
    run(croo_cmd; wait=true)
    @info string("running ", croo_cmd, " done!")
    cp(rna_metadata_json_path, joinpath(output_dir, basename(rna_metadata_json_path)))
end
@info "collecting outputs done!"

# extract RSEM quantification results
@info "start extracting RSEM quantification results ..."
genes_pattern = r"_anno_rsem\.genes\.results$"
isoforms_pattern = r"_anno_rsem\.isoforms\.results$"
group_pattern = r"^rep[0-9]+|_anno_rsem\.(genes|isoforms)\.results$"
replicate_pattern = r"(?<replicate>^rep[0-9]+)"

genes_output_dir = joinpath(args["output_dir"][1], "genes")
isoforms_output_dir = joinpath(args["output_dir"][1], "isoforms")

mkpath(genes_output_dir)
mkpath(isoforms_output_dir)

genes_files = YRUtils.BaseUtils.list_files(args["output_dir"][1], genes_pattern; recursive=true, full_name=true)
isoforms_files = YRUtils.BaseUtils.list_files(args["output_dir"][1], isoforms_pattern; recursive=true, full_name=true)

genes_df = DataFrame(old_file=genes_files,
    old_basename=basename.(genes_files),
    group=replace.(basename.(genes_files), group_pattern => ""),
    replicate=[m["replicate"] for m in match.(replicate_pattern, basename.(genes_files))])
genes_df.sample = string.(genes_df.group, "_", genes_df.replicate)
genes_df.new_basename = string.(genes_df.sample, "_anno_rsem.genes.results")
genes_df.new_file = joinpath.(genes_output_dir, genes_df.new_basename)
cp.(genes_df.old_file, genes_df.new_file)
CSV.write(joinpath(genes_output_dir, "genes.sample_sheet.tsv"), genes_df; delim="\t", writeheader=true)

isoforms_df = DataFrame(old_file=isoforms_files,
    old_basename=basename.(isoforms_files),
    group=replace.(basename.(isoforms_files), group_pattern => ""),
    replicate=[m["replicate"] for m in match.(replicate_pattern, basename.(isoforms_files))])
isoforms_df.sample = string.(isoforms_df.group, "_", isoforms_df.replicate)
isoforms_df.new_basename = string.(isoforms_df.sample, "_anno_rsem.isoforms.results")
isoforms_df.new_file = joinpath.(isoforms_output_dir, isoforms_df.new_basename)
cp.(isoforms_df.old_file, isoforms_df.new_file)
CSV.write(joinpath(isoforms_output_dir, "isoforms.sample_sheet.tsv"), isoforms_df; delim="\t", writeheader=true)
@info "extracting RSEM quantification results done!"

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
