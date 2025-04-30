
echo "## add environment variable rsRNASP_RNA_HOME into ~/.bashrc"
echo "## please restart the terminal"
rsRNASP_RNA_HOME=`pwd`
user_home=`echo ~`
envsetnot=`env | awk -F= '{print $1}' | grep "rsRNASP_RNA_HOME"`
if [ -z "$envsetnot" ];then
        cat   << EOF >> ${user_home}/.bashrc
export rsRNASP_RNA_HOME=${rsRNASP_RNA_HOME}
EOF

fi


