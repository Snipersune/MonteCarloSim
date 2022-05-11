#!/bin/bash

ls data/064_* > data/064all_files
ls data/256_* > data/256all_files
ls data/1024_* > data/1024all_files

VAR64=$(echo $(cat data/064all_files))
VAR256=$(echo $(cat data/256all_files))
VAR1024=$(echo $(cat data/1024all_files))

tail -q -n 1 ${VAR64} > data/_064_data
tail -q -n 1 ${VAR256} > data/_256_data
tail -q -n 1 ${VAR1024} > data/_1024_data

rm data/064all_files
rm data/256all_files
rm data/1024all_files
