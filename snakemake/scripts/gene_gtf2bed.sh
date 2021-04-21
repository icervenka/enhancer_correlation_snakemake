#!/bin/bash

INPUT=$1

awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$16,$7}}' | tr -d '";' $1 
