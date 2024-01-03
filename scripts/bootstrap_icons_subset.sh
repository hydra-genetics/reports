#!/bin/bash

# Extracts glyphs from a font using fonttools.
# fontpath should be the path to the font file and unicodes should be a
# comma-separated list of glyphs in unicode representation, or unicode
# ranges. See the fonttools documentation for more information.
#
# A base64 encoded version of the subset font file will be printed on
# stdout, and various useful information will be printed to stderr.

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 unicodes"
    exit 1
fi

FONTURL="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.2/font/fonts/bootstrap-icons.woff2"
FONTFILE=$(basename "${FONTURL}")
EXT="${FONTFILE##*.}"

FONTPATH="${TMPDIR:-/tmp}/bootstrap-icons.${EXT}"
OUTPATH="${TMPDIR:-/tmp}/${FONTFILE}.subset.${EXT}"

UNICODES="$1"

curl --silent --output "${FONTPATH}" "${FONTURL}"

fonttools subset "${FONTPATH}" \
    --glyph-names \
    --flavor="${EXT}" \
    --unicodes="${UNICODES}" \
    --output-file="${OUTPATH}"
fonttools ttx -t cmap --flavor="${EXT}" -o "${OUTPATH}.xml" "${OUTPATH}"

cat "${OUTPATH}.xml" >&2
echo >&2
base64 -w0 "${OUTPATH}"
echo -e "\n"
file "${OUTPATH}" >&2

rm "${FONTPATH}"
rm "${OUTPATH}"
rm "${OUTPATH}.xml"
