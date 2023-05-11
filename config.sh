while test $# -gt 0
do
    case "$1" in
	ceci)
            source load_ceci_modules.sh
            ln -s $(pwd) ~/classicalDFT
            ;;
        *) echo "bad option $1"
            ;;
    esac
    shift
done

#! /bin/sh
root_dir=$(pwd)
parent_dir="$(dirname "$root_dir")"

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release $1 ..
cd $root_dir

mkdir -p debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug $1 ..
cd $root_dir
