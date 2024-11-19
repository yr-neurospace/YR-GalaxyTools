using CSV, DataFrames

TARGET_JOBS = ("encode_bulk_rna_seq_pipeline", "encode_bulk_atac_seq_pipeline", "encode_bulk_chip_seq_pipeline")
TSV_HEADERS = ["tool_id", "tool_version", "destination_id", "handler", "state", "job_runner_name", "sum", "user_id"]
WELCOME_PAGE_PATH = "/srv/galaxy/var/config/welcome.html"

stdout_io_buffer = IOBuffer()
run(pipeline(`sudo --stdin --user galaxy gxadmin tsvquery queue-overview`;
        stdin=IOBuffer("rootpw159*"), stdout=stdout_io_buffer); wait=true)
df = CSV.read(IOBuffer(String(take!(stdout_io_buffer))), DataFrame; delim="\t", header=TSV_HEADERS, types=String)
subset!(df, :tool_id => x -> in.(x, Ref(TARGET_JOBS)))

raw_rows = [[r[i] for i in propertynames(df)] for r in eachrow(df)]
td_rows = [string.("<td>", r, "</td>") for r in raw_rows]
tr_rows = [string("<tr>", join(r, ""), "</tr>") for r in td_rows]

old_html_str = read(WELCOME_PAGE_PATH, String)
if !isnothing(match(Regex("<tbody>.*</tbody>", "s"), old_html_str))
    new_html_str = replace(old_html_str, Regex("<tbody>.*</tbody>", "s") => string("<tbody>", join(tr_rows, ""), "</tbody>"))
    write(WELCOME_PAGE_PATH, new_html_str)
else
    @error string("cannot find <tbody>...</tbody> in ", WELCOME_PAGE_PATH)
end