BEGIN {
}

/background in box/ {
    i = 0;
}

/signal scaling:/ {
    scale[i] = $3;
}

/CL\(s\):/ {
    cls[i] = $2;
    printf("%f %f \n",scale[i],cls[i]);
    i++;
}

END {
}

