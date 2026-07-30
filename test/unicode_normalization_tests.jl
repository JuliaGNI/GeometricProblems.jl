using Test
using Unicode: normalize

# Every source file in this repository is NFC-normalized.
#
# Accented identifiers have two legal spellings: `ū` as the single codepoint U+016B, or `u` followed
# by U+0304 COMBINING MACRON. Julia's parser NFC-normalizes identifiers, so the two are the *same
# binding* — a file mixing them compiles and runs identically, and no test of behaviour can tell them
# apart. What they are not is the same *text*: `grep ū` typed one way silently finds nothing while the
# identifier sits plainly on screen, and an exact-string edit fails for no visible reason. Files drift
# between the forms without anyone touching a character, because macOS filesystem APIs hand back
# decomposed text and some editors recompose on save.
#
# This asserts the normal form of the file rather than any particular character, so there is no list
# of accented identifiers to keep up to date. Characters with no precomposed form pass untouched:
# `q̇` is `q` followed by U+0307 and there is no single codepoint for it, so it is already in normal
# form. `obsolete/` is excluded — it is kept for reference and not edited.

const REPO_ROOT = dirname(@__DIR__)
const SKIP_DIRS = ("obsolete", "build", "node_modules")
const EXTENSIONS = (".jl", ".md", ".toml", ".yml", ".bib")

"Every text file in the repository whose normal form this test is responsible for."
function source_files(root = REPO_ROOT)
    files = String[]
    for (dir, subdirs, names) in walkdir(root)
        # prune rather than filter afterwards, so that `.git` and friends are never descended into
        filter!(d -> d ∉ SKIP_DIRS && !startswith(d, "."), subdirs)
        for name in names
            any(ext -> endswith(name, ext), EXTENSIONS) && push!(files, joinpath(dir, name))
        end
    end
    return sort!(files)
end

@testset "$(rpad("Unicode normalization (NFC)", 80))" begin
    files = source_files()

    # A walk that found nothing would make the assertion below pass vacuously.
    @test !isempty(files)

    offenders = String[]
    for file in files
        text = read(file, String)
        normalize(text, :NFC) == text || push!(offenders, relpath(file, REPO_ROOT))
    end

    # Named rather than counted, so that a failure says which file to run `normalize` over.
    @test offenders == String[]
end
