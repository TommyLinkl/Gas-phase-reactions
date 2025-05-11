#!/bin/bash

# Loop through all files ending with _ini_fin.xyz in current directory
for file in *_ini_fin.xyz; do
    # Check that the file exists (avoid error if no match)
    [ -e "$file" ] || continue

    # Remove extension and append interpolation suffix
    base="${file%.xyz}"
    output="${base}_interpolation.xyz"

    # Run geodesic_interpolate with 20 images
    geodesic_interpolate "$file" --output "$output" --nimages 20
done

