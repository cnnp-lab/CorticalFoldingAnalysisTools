function lobe_label_ds=label2lobe(label_ds,LUT_lobes)
lobe_label_ds=zeros(size(label_ds));
ul=unique(label_ds);
for ui=2:length(ul)%starts frm 2 as 1 is the the 0 label
    lid=find(LUT_lobes(:,1)==ul(ui));
    lobe_label_ds(label_ds==ul(ui))=LUT_lobes(lid,2);
end