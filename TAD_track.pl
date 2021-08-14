#!usr/bin/perl
$file=$ARGV[0];
$out=$ARGV[1];
open FILE,"$file";
@name=split(/\//,$file);
$length=scalar(@name)-1;
$file=~s/$name[$length]/$out/g;
open OUT,">$file.temp";
$i=1;
$line=<FILE>;
while($line=<FILE>){
    chomp $line;
    @info=split(/\t/,$line);
    $print1=$info[1]+1000;
    $print2=$info[5]-1000;
    $print3=$info[2]-$info[1];
    $print4=$info[5]-$info[4];
    $print5=$info[4]-$info[1];
    print OUT "$info[0]\t$info[1]\t$info[2]\n";
    $i++;
}
system("sort -k1,1 -k2,2n $file.temp > $file");
system("rm $file.temp");
#system("sed "1i\track name=$out description=$out type=bigBed" $file");
$bb=$file;
$bb=~s/bed/bb/g;
system("~/bedToBigBed $file /home/lingshi/genome/hg38_p2/hg38_2.chrom.sizes $bb");
system("rm $file");
