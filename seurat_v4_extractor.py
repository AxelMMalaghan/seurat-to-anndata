from seurat_extractor import SeuratExtractor


class SeuratV4Extractor(SeuratExtractor):

    def _extract_assay(self, data, assay_name, layer_name):

        assay_pointer = data['assays']['data'][assay_name]['data']

        X = assay_pointer[layer_name].T
        var_names = assay_pointer['features']['data'] if 'features' in assay_pointer else None
        return var_names, X