function output = gmb_tbl_nanblank_V0(values)
    mask = isnan(values);
    if nnz(mask)
      output = string(values);
      output(mask) = "";
      output = char(output);
    else
        output = values;
    end
end