SHELL=/bin/bash

REMOTE_BB=https://resources.altius.org/~areynolds/public/masterlist_DHSs_733samples_WM20180608_all_mean_signal_colorsMax.bed.unc.b12.lf.AR20210927.withBiosampleCounts.bb
LOCAL_PREFIX=masterlist_DHSs_733samples_WM20180608_all_mean_signal_colorsMax.bed.unc.b12.lf.AR20210927.withBiosampleCounts
LOCAL_BB=${LOCAL_PREFIX}.bb
LOCAL_BED12=${LOCAL_PREFIX}.bed

REMOTE_BB_TO_BED_OSX=https://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/bigBedToBed
LOCAL_BB_TO_BED_OSX=bigBedToBed

BED12_TO_SAM=bed12_to_sam.py

all: prep extract split convert check package upload

prep:
	wget -qO- ${REMOTE_BB} > ${LOCAL_BB}
	wget -qO- ${REMOTE_BB_TO_BED_OSX} > ${LOCAL_BB_TO_BED_OSX}
	chmod +x ${LOCAL_BB_TO_BED_OSX}

clean:
	rm -f ${LOCAL_BB} ${LOCAL_BB_TO_BED_OSX} ${LOCAL_BED12} ${LOCAL_PREFIX}.chr*.bed ${LOCAL_PREFIX}.chr*.sam ${LOCAL_PREFIX}.chr*.bam ${LOCAL_PREFIX}.chr*.bam.bai

extract:
	${PWD}/${LOCAL_BB_TO_BED_OSX} ${LOCAL_BB} ${LOCAL_BED12}

split:
	for i in `bedextract --list-chr ${LOCAL_BED12}`; do echo $${i}; bedextract $${i} ${LOCAL_BED12} > ${LOCAL_PREFIX}.$${i}.bed; done

convert:
	for i in `bedextract --list-chr ${LOCAL_BED12}`; do echo $${i}; ${PWD}/${BED12_TO_SAM} ${LOCAL_PREFIX}.$${i}.bed ${LOCAL_PREFIX}.$${i}.sam; done

check:
	for i in `bedextract --list-chr ${LOCAL_BED12}`; do echo $${i}; samtools quickcheck -v ${LOCAL_PREFIX}.$${i}.sam; done

package:
	for i in `bedextract --list-chr ${LOCAL_BED12}`; do echo $${i}; samtools view -bS ${LOCAL_PREFIX}.$${i}.sam > ${LOCAL_PREFIX}.$${i}.bam; samtools index ${LOCAL_PREFIX}.$${i}.bam; done

upload:
	aws s3 cp . s3://areynolds-us-west-2 --recursive --exclude "*" --include "${LOCAL_PREFIX}.chr*.bam"
	aws s3 cp . s3://areynolds-us-west-2 --recursive --exclude "*" --include "${LOCAL_PREFIX}.chr*.bam.bai"