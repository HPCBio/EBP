# Software

The software and documentation files are provided by the vendor.
https://www.universalsequencing.com/analysis-tools

You should receive these files from the download link that the vendors provide.

Download them locally to your laptop.

<pre>
TELL Seq Software User Guide (Tell-Link) (1).pdf  
TELL Seq Software User Guide (Tell-Sort).pdf
TELL-Seq Software User Guide (Tell-Read).pdf      
conversion_tool (4).tar
tellink-v1.0.2.tar
tellread-v1.0.2.tar
tellsort-v1.0.2.tar
</pre>

# Dependencies

The actual analysis tools are provided as docker containers. Therefore, there are no pre-requisites. 

Make sure you can install and execute docker container using Docker or Singularity

In our case, we need Singularity already installed and configured on the cluster

- singularity version 3.4.1


# Installation 

This step should be straightforward as containers are by definition "out-of-the-box" toolkits.

However, the path to the installation folder is determined by how Docker or Singularity was setup on the cluster.

These were the steps for installing the software so that it could work with Singularity on the biotech server

## Untar vendor files

<pre>

$ scp * grendon@compute-4.biotech.illinois.edu:/home/a-m/grendon/TELL-Seq/v1.0.2/

$ ssh grendon@compute-4.biotech.illinois.edu

$ cd /home/a-m/grendon/TELL-Seq/v1.0.2/

$ ls 
TELL Seq Software User Guide (Tell-Link) (1).pdf  
TELL Seq Software User Guide (Tell-Sort).pdf
TELL-Seq Software User Guide (Tell-Read).pdf      
conversion_tool (4).tar
tellink-v1.0.2.tar
tellread-v1.0.2.tar
tellsort-v1.0.2.tar

$ tar -xvf tellread-v1.0.2.tar
tellread-release/
tellread-release/generateGenomeIndexBed.sh
tellread-release/run_tellread_fq.sh
tellread-release/docker-tellread
tellread-release/run_tellread.sh


$ tar -xvf tellink-v1.0.2.tar
tellink-release/
tellink-release/docker-tellink
tellink-release/run_tellink.sh

</pre>

## Install containers with Singularity

<pre>

$ module load singularity/3.4.1
$ cd tellread-release/
$ singularity build tellread.simg docker-archive://docker-tellread
 

INFO:    Starting build...
Getting image source signatures
Copying blob sha256:c8be1b8f4d60d99c281fc2db75e0f56df42a83ad2f0b091621ce19357e19d853
 62.54 MiB / 62.54 MiB [====================================================] 7s
Copying blob sha256:977183d4e9995d9cd5ffdfc0f29e911ec9de777bcb0f507895daa1068477f76f
 968.00 KiB / 968.00 KiB [==================================================] 0s
Copying blob sha256:6597da2e2e52f4d438ad49a14ca79324f130a9ea08745505aa174a8db51cb79d
 15.50 KiB / 15.50 KiB [====================================================] 0s
Copying blob sha256:16542a8fc3be1bfaff6ed1daa7922e7c3b47b6c3a8d98b7fca58b9517bb99b75
 3.00 KiB / 3.00 KiB [======================================================] 0s
Copying blob sha256:5866d7158bc4687bfa5cb75159af6aaef9937abaa3037d1ec83042cde8789e31
 2.83 GiB / 2.83 GiB [===================================================] 3m54s
Copying blob sha256:1d313c5679ea2142dec7fa9a73d2d1fb7fc50e8ce89d0493e1aefeec3b877eec
 18.92 MiB / 18.92 MiB [====================================================] 1s
Copying config sha256:277d5ea02452abccb906234ea3a971b43727a3dc0b5617e3c3a9d0b231c58de1
 4.44 KiB / 4.44 KiB [======================================================] 0s
Writing manifest to image destination
Storing signatures
2020/10/05 13:10:08  info unpack layer: sha256:5bed26d33875e6da1d9ff9a1054c5fef3bbeb22ee979e14b72acf72528de007b
2020/10/05 13:10:08  warn xatt{etc/gshadow} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:08  warn xatt{etc/shadow} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:08  warn xatt{run/utmp} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:08  warn xatt{sbin/pam_extrausers_chkpwd} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:08  warn xatt{sbin/unix_chkpwd} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:08  warn xatt{usr/bin/chage} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:08  warn xatt{usr/bin/expiry} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:08  warn xatt{usr/bin/wall} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:09  warn xatt{var/cache/apt/archives/partial} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:09  warn xatt{var/local} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:09  warn xatt{var/log/apt/term.log} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:09  warn xatt{var/log/btmp} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:09  warn xatt{var/log/lastlog} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:09  warn xatt{var/log/wtmp} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:09  warn xatt{var/mail} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:09  info unpack layer: sha256:f11b29a9c7306674a9479158c1b4259938af11b97359d9ac02030cc1095e9ed1
2020/10/05 13:10:09  info unpack layer: sha256:930bda195c84cf132344bf38edcad255317382f910503fef234a9ce3bff0f4dd
2020/10/05 13:10:09  info unpack layer: sha256:78bf9a5ad49e4ae42a83f4995ade4efc096f78fd38299cf05bc041e8cdda2a36
2020/10/05 13:10:09  info unpack layer: sha256:31e1e89520eedcc1070c9661091a5916006d88a287dae8cd0dcf49ac603c5dfb
2020/10/05 13:10:52  info unpack layer: sha256:b26524500544222c6113cda3f3c31520fb7bc053819ea50043aabac860b5fd3b
2020/10/05 13:10:52  warn xatt{var/cache/apt/archives/partial} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:10:52  warn xatt{var/log/apt/term.log} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
INFO:    Creating SIF file...
INFO:    Build complete: tellread.simg


$ cd ../tellink-release/

$ singularity build tellink.simg docker-archive://docker-tellink

INFO:    Starting build...
Getting image source signatures
Copying blob sha256:c8be1b8f4d60d99c281fc2db75e0f56df42a83ad2f0b091621ce19357e19d853
 62.54 MiB / 62.54 MiB [====================================================] 7s
Copying blob sha256:977183d4e9995d9cd5ffdfc0f29e911ec9de777bcb0f507895daa1068477f76f
 968.00 KiB / 968.00 KiB [==================================================] 0s
Copying blob sha256:6597da2e2e52f4d438ad49a14ca79324f130a9ea08745505aa174a8db51cb79d
 15.50 KiB / 15.50 KiB [====================================================] 0s
Copying blob sha256:16542a8fc3be1bfaff6ed1daa7922e7c3b47b6c3a8d98b7fca58b9517bb99b75
 3.00 KiB / 3.00 KiB [======================================================] 0s
Copying blob sha256:575e42fb2af20e621101ea32ce13296ada8788cfc8c61aaf27cab1082d9bb5e7
 981.87 MiB / 981.87 MiB [===============================================] 1m31s
Copying blob sha256:5a42686574db8eb8f3977decd692d2af480a70c103ff028108264d9b6fcc1e26
 18.92 MiB / 18.92 MiB [====================================================] 1s
Copying config sha256:859c13da8fe36a0dee5498d9283bd6c2b633ffd71ee8fe8777d8e788c00d33c9
 4.48 KiB / 4.48 KiB [======================================================] 0s
Writing manifest to image destination
Storing signatures
2020/10/05 13:19:04  info unpack layer: sha256:5bed26d33875e6da1d9ff9a1054c5fef3bbeb22ee979e14b72acf72528de007b
2020/10/05 13:19:05  warn xatt{etc/gshadow} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:05  warn xatt{etc/shadow} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:05  warn xatt{run/utmp} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:05  warn xatt{sbin/pam_extrausers_chkpwd} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:05  warn xatt{sbin/unix_chkpwd} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:05  warn xatt{usr/bin/chage} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:05  warn xatt{usr/bin/expiry} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:05  warn xatt{usr/bin/wall} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:06  warn xatt{var/cache/apt/archives/partial} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:06  warn xatt{var/local} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:06  warn xatt{var/log/apt/term.log} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:06  warn xatt{var/log/btmp} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:06  warn xatt{var/log/lastlog} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:06  warn xatt{var/log/wtmp} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:06  warn xatt{var/mail} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:06  info unpack layer: sha256:f11b29a9c7306674a9479158c1b4259938af11b97359d9ac02030cc1095e9ed1
2020/10/05 13:19:06  info unpack layer: sha256:930bda195c84cf132344bf38edcad255317382f910503fef234a9ce3bff0f4dd
2020/10/05 13:19:06  info unpack layer: sha256:78bf9a5ad49e4ae42a83f4995ade4efc096f78fd38299cf05bc041e8cdda2a36
2020/10/05 13:19:06  info unpack layer: sha256:fce278f53e94a357096edb60d5ceb325aa8cab04da81cf9a0c3e038a57d87f77
2020/10/05 13:19:20  info unpack layer: sha256:5d113128ceecf9c5e0f03c7a31451a05ad1701f17289ac70601449154a7b553f
2020/10/05 13:19:21  warn xatt{var/cache/apt/archives/partial} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
2020/10/05 13:19:21  warn xatt{var/log/apt/term.log} ignoring ENOTSUP on setxattr "user.rootlesscontainers"
INFO:    Creating SIF file...
INFO:    Build complete: tellink.simg

</pre>

## The folder structure

After installing the tools using Singularity, the folder should look like this:


<pre>

$ cd /home/a-m/grendon/TELL-Seq/v1.0.2/

$ tree
.
|-- TELL\ Seq\ Software\ User\ Guide\ (Tell-Link)\ (1).pdf
|-- TELL\ Seq\ Software\ User\ Guide\ (Tell-Sort).pdf
|-- TELL-Seq\ Software\ User\ Guide\ (Tell-Read).pdf
|-- conversion_tool\ (4).tar
|-- tellink-release
|   |-- docker-tellink
|   |-- run_tellink.sh
|   `-- tellink.simg
|-- tellink-v1.0.2.tar
|-- tellread-release
|   |-- docker-tellread
|   |-- generateGenomeIndexBed.sh
|   |-- run_tellread.sh
|   |-- run_tellread_fq.sh
|   `-- tellread.simg
|-- tellread-v1.0.2.tar
`-- tellsort-v1.0.2.tar

</pre>
