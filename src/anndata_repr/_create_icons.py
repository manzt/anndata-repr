import math
import pathlib
import uuid

def get_display(adata):

    # SIZING
    width = 1000
    height = 1000

    n_rows_x = adata.X.shape[0]
    n_cols_x = adata.X.shape[1]
    n_rows_var = adata.var.shape[1]
    n_cols_obs = adata.obs.shape[1]

    n_layers_x = len(adata.layers)
    n_layers_obsm = len(adata.obsm)
    n_layers_obsp = len(adata.obsp)
    n_layers_varm = len(adata.varm)
    n_layers_varp = len(adata.varp)

    margin_side = 20
    margin_layer = 2
    margin_blocks = 10

    # Get the size of x
    size_rows_x, size_cols_x = get_size_x(n_rows_x, n_cols_x, width, height)

    # Set default sizing for width/height of additional blocks
    size_cols_obs = 20
    size_rows_var = 20
    size_cols_obsm = 10
    size_rows_varm = 10

    # X/Y offsets
    x_base_obsp = margin_side
    x_base_obsm = x_base_obsp + size_rows_x + margin_blocks if n_layers_obsp > 0 else x_base_obsp
    x_base = x_base_obsm + size_cols_obsm + margin_blocks if n_layers_obsm > 0 else x_base_obsm
    x_base_obs = x_base + size_cols_x + margin_blocks

    x_total = x_base_obs + size_cols_obs + margin_side

    y_base_var = margin_side
    y_base = y_base_var + size_rows_var + margin_blocks
    y_base_varm = y_base + size_rows_x + margin_blocks
    y_base_varp = y_base_varm + size_rows_varm + margin_blocks if n_layers_varm > 0 else y_base_varm

    y_extra = 100 if len(adata.obs) > 0 else 0

    y_total = y_base_varp + size_cols_x + margin_side + y_extra if n_layers_varp > 0 else y_base_varp + margin_side + y_extra


    # SVG

    # TOP LEVEL
    style = (
        '''<style>
            .cls-2,.cls-21{fill:#efc41c;}
            .cls-11,.cls-2,.cls-4,.cls-5,.cls-6,.cls-7,.cls-8,.cls-9{stroke:#010101;stroke-miterlimit:10;stroke-width:0.5px;}
            .cls-3{fill:#fff;}
            .cls-4{fill:#f15c5a;}
            .cls-5{fill:#965ba5;}
            .cls-6{fill:#194c61;}
            .cls-7{fill:#ef9021;}
            .cls-10,.cls-8{fill:#3c8b53;}
            .cls-9{fill:#4fba6f;}
            .cls-11,.cls-25{fill:#2c96c0;}
            .uns-braces{fill:#010101}
            .uns-name{fill:#231f20}
            .uns-stroke{stroke:white;stroke-width:2px}
            .var-line{fill:#2c96c0;}
            .obs-line{fill:#efc41c;}
            .cls-17{fill:#4899d4;}
            #var-names-horizontal .var-line{mask:url(#mask-var-names-horizontal)}
            #obs-names-horizontal .obs-line{mask:url(#mask-obs-names-horizontal)}
            #var-names-vertical .var-line{mask:url(#mask-var-names-vertical)}
            #obs-names-vertical .obs-line{mask:url(#mask-obs-names-vertical)}
        </style>'''
    )

    uns_specification = (
        '''
        <g id="uns" class=block">
            <path id="uns-braces-left" d="M24.58916,159.27075q0-.31786.02-.84668c.01319-.35254.02637-.7207.04-1.10644.0127-.38574.02637-.76367.03955-1.13672.01319-.37109.02-.69043.02-.95606a7.28108,7.28108,0,0,0-.15967-1.6748,2.78178,2.78178,0,0,0-.45849-1.04687,1.608,1.608,0,0,0-.74707-.54786,2.99014,2.99014,0,0,0-1.00684-.15918H20.62237v-3.62793h1.67431a2.26225,2.26225,0,0,0,1.874-.76757,4.01,4.01,0,0,0,.61767-2.502q0-.23877-.00976-.708c-.00684-.31153-.01661-.6377-.02979-.97657q-.02051-.50829-.03027-.94726-.01026-.438-.00977-.59766a7.23611,7.23611,0,0,1,1.95362-5.50195q1.95336-1.853,5.84082-1.85352h3.30908v3.6875H31.80547a2.29678,2.29678,0,0,0-1.94385.81738,3.99329,3.99329,0,0,0-.62793,2.45215q0,.2798.01026.72754.00952.44825.02978.92676.01978.479.02979.89746c.00683.27832.01025.48535.01025.61719a6.66862,6.66862,0,0,1-.76758,3.22949,6.21961,6.21961,0,0,1-2.123,2.293,5.262,5.262,0,0,1,1.12646.98633,6.5702,6.5702,0,0,1,.88672,1.33594,7.32034,7.32034,0,0,1,.58838,1.63476,8.1971,8.1971,0,0,1,.209,1.88379c0,.14648-.00341.40137-.00976.7666q-.01026.54932-.02979,1.11719-.02051.56689-.03027,1.03613c-.00684.3125-.00977.49512-.00977.54785a4.30526,4.30526,0,0,0,.877,3.02051,3.84045,3.84045,0,0,0,2.89062.92676h2.85059v3.6084H33.47979q-.41821,0-1.24561-.04981-.82764-.0498-1.77441-.19922a14.75574,14.75574,0,0,1-1.86377-.41894,4.604,4.604,0,0,1-1.51514-.72754,5.73374,5.73374,0,0,1-1.87353-2.542A9.90217,9.90217,0,0,1,24.58916,159.27075Z"/>
            <path id="uns-braces-dots" d="M39.4002,158.1145a2.823,2.823,0,0,1,.20947-1.08593,2.514,2.514,0,0,1,.59815-.877,2.9987,2.9987,0,0,1,.90673-.58789,2.93511,2.93511,0,0,1,1.15625-.21972,3.03122,3.03122,0,0,1,1.15625.21972,2.82393,2.82393,0,0,1,.917.58789,2.77809,2.77809,0,0,1,.59814,3.04981,2.82037,2.82037,0,0,1-.59814.88672,2.74,2.74,0,0,1-.917.59863,3.04921,3.04921,0,0,1-1.15625.21875,2.95229,2.95229,0,0,1-1.15625-.21875,2.90362,2.90362,0,0,1-.90673-.59863,2.60836,2.60836,0,0,1-.59815-.88672A2.82435,2.82435,0,0,1,39.4002,158.1145Zm6.77783,0a2.83179,2.83179,0,0,1,.209-1.08593,2.514,2.514,0,0,1,.59815-.877,3.00716,3.00716,0,0,1,.90723-.58789,2.93161,2.93161,0,0,1,1.15625-.21972,3.027,3.027,0,0,1,1.15576.21972,2.82389,2.82389,0,0,1,.917.58789,2.77809,2.77809,0,0,1,.59814,3.04981,2.82037,2.82037,0,0,1-.59814.88672,2.74,2.74,0,0,1-.917.59863,3.04493,3.04493,0,0,1-1.15576.21875,2.94877,2.94877,0,0,1-1.15625-.21875,2.91157,2.91157,0,0,1-.90723-.59863,2.60836,2.60836,0,0,1-.59815-.88672A2.8332,2.8332,0,0,1,46.178,158.1145Zm6.73779,0a2.83179,2.83179,0,0,1,.209-1.08593,2.57758,2.57758,0,0,1,.58838-.877,2.84483,2.84483,0,0,1,.90673-.58789,2.95539,2.95539,0,0,1,1.14649-.21972,3.03126,3.03126,0,0,1,1.15625.21972,2.95086,2.95086,0,0,1,.92676.58789,2.6349,2.6349,0,0,1,.60791.877,2.75832,2.75832,0,0,1-.60791,3.05957,2.85963,2.85963,0,0,1-.92676.59863,3.04925,3.04925,0,0,1-1.15625.21875,2.97283,2.97283,0,0,1-1.14649-.21875,2.75929,2.75929,0,0,1-.90673-.59863,2.67679,2.67679,0,0,1-.58838-.88672A2.8332,2.8332,0,0,1,52.91582,158.1145Z"/>
            <path id="uns-braces-right" d="M74.26494,159.27075a9.91264,9.91264,0,0,1-.61767,3.61817,5.73463,5.73463,0,0,1-1.874,2.542,4.6011,4.6011,0,0,1-1.51464.72754,14.78235,14.78235,0,0,1-1.86426.41894q-.947.14941-1.77393.19922-.82763.04981-1.24609.04981h-2.292V163.218h2.85058a3.84,3.84,0,0,0,2.89014-.92676,4.30521,4.30521,0,0,0,.87744-3.02051c0-.05273-.00391-.23535-.01025-.54785q-.01026-.46875-.02979-1.03613-.0205-.56836-.02978-1.11719c-.00684-.36523-.01026-.62012-.01026-.7666a8.16977,8.16977,0,0,1,.20948-1.88379,7.63525,7.63525,0,0,1,.57812-1.63476,6.034,6.034,0,0,1,.88721-1.33594,5.454,5.454,0,0,1,1.13623-.98633,6.21159,6.21159,0,0,1-2.123-2.293,6.66312,6.66312,0,0,1-.76758-3.22949c0-.13184.00293-.33887.00977-.61719q.01024-.419.03027-.89746.019-.479.02979-.92676.00951-.44825.00976-.72754a3.98994,3.98994,0,0,0-.62793-2.45215,2.29516,2.29516,0,0,0-1.94336-.81738H63.04229v-3.6875h3.30908q3.887,0,5.84082,1.85352a7.238,7.238,0,0,1,1.95361,5.50195q0,.15968-.01025.59766-.01026.43945-.02979.94726c-.01367.33887-.02343.665-.02978.97657q-.01026.46875-.01025.708a4.00652,4.00652,0,0,0,.61816,2.502,2.26121,2.26121,0,0,0,1.87353.76757h1.67481v3.62793H76.51787a2.963,2.963,0,0,0-1.0166.15918,1.61537,1.61537,0,0,0-.73779.54786,2.78952,2.78952,0,0,0-.4585,1.04687,7.28493,7.28493,0,0,0-.15918,1.6748c0,.26563.00635.585.01953.95606.01319.373.02686.751.04,1.13672s.02637.7539.04,1.10644Q74.26446,158.95288,74.26494,159.27075Z"/>
            <path id="uns-name-u" d="M46.89336,134.05884l.0083,3.08008a1.07024,1.07024,0,0,0,.02149.20117q.01977.10547.0498.21973l-.813.00195a.1628.1628,0,0,1-.05908-.09863.54772.54772,0,0,1-.01465-.11719l-.00049-.11328a.94536.94536,0,0,1-.186.20215,1.288,1.288,0,0,1-.23779.15039,1.33707,1.33707,0,0,1-.26123.0957,1.05915,1.05915,0,0,1-.2583.03418,1.40279,1.40279,0,0,1-.49122-.08008.86672.86672,0,0,1-.36669-.25488,1.17934,1.17934,0,0,1-.23047-.45313,2.51078,2.51078,0,0,1-.08057-.66992l-.00586-2.19043.82764-.002.00586,2.18945a.99945.99945,0,0,0,.13916.58789.42983.42983,0,0,0,.36474.1836.68125.68125,0,0,0,.29444-.06446.71848.71848,0,0,0,.23779-.17871.84461.84461,0,0,0,.15869-.27344,1.00608,1.00608,0,0,0,.05762-.3496l-.00586-2.09864Z"/>
            <path id="uns-name-n" d="M47.79668,137.55689l-.00976-3.501.84619-.002.001.35839a1.44724,1.44724,0,0,1,.42383-.31738,1.10427,1.10427,0,0,1,.49024-.12109,1.36743,1.36743,0,0,1,.48.08008.88354.88354,0,0,1,.36865.25293,1.1782,1.1782,0,0,1,.23584.45019,2.35643,2.35643,0,0,1,.08447.668l.00586,2.124-.82763.002-.00586-2.11621a1.47975,1.47975,0,0,0-.0376-.34864.72162.72162,0,0,0-.10352-.23925.41475.41475,0,0,0-.15771-.13672.45765.45765,0,0,0-.20557-.043.7069.7069,0,0,0-.29394.0625.69985.69985,0,0,0-.2378.17676.85516.85516,0,0,0-.15869.27343,1.018,1.018,0,0,0-.05762.34961l.00586,2.02539Z"/>
            <path id="uns-name-s" d="M54.01543,135.177a.06119.06119,0,0,1-.04541-.01562.12305.12305,0,0,1-.02588-.04.53849.53849,0,0,1-.0166-.05176.27425.27425,0,0,0-.01855-.04883.95758.95758,0,0,0-.165-.15527,1.16635,1.16635,0,0,0-.21826-.12988,1.40357,1.40357,0,0,0-.25293-.08692,1.1324,1.1324,0,0,0-.26563-.03222,1.38254,1.38254,0,0,0-.21582.0166.7157.7157,0,0,0-.18164.05371.339.339,0,0,0-.124.0957.22409.22409,0,0,0-.04492.14356.249.249,0,0,0,.0498.15918.46265.46265,0,0,0,.148.11621,1.586,1.586,0,0,0,.24951.10058c.1001.03125.21875.06446.356.10157a2.22561,2.22561,0,0,1,.92187.44726.84089.84089,0,0,1,.29688.62891.98842.98842,0,0,1-.10547.45117,1.18346,1.18346,0,0,1-.29834.376,1.4596,1.4596,0,0,1-.47022.25976,1.88465,1.88465,0,0,1-.61669.09668,2.316,2.316,0,0,1-.81885-.13867,2.10663,2.10663,0,0,1-.70655-.4502l.4336-.78906a.08722.08722,0,0,1,.0664.05371c.01465.03028.03272.07227.05518.126a.74224.74224,0,0,0,.169.18457,1.53352,1.53352,0,0,0,.24609.16211,1.45508,1.45508,0,0,0,.28369.11621.96614.96614,0,0,0,.28418.043,1.01266,1.01266,0,0,0,.46436-.08692.2768.2768,0,0,0,.16455-.25878.31519.31519,0,0,0-.03711-.15332.40012.40012,0,0,0-.11573-.123.97049.97049,0,0,0-.20751-.10547c-.08545-.03223-.19043-.06738-.315-.10352-.06836-.01855-.13623-.03808-.20361-.05761s-.13525-.043-.20312-.06934a2.16868,2.16868,0,0,1-.3667-.17285,1.20388,1.20388,0,0,1-.27344-.21777.83353.83353,0,0,1-.17139-.27149.94521.94521,0,0,1-.05957-.33887.76563.76563,0,0,1,.10352-.374,1.16939,1.16939,0,0,1,.28613-.33594,1.52974,1.52974,0,0,1,.4375-.24218,1.59669,1.59669,0,0,1,.55762-.09571,1.91722,1.91722,0,0,1,.74219.14063,1.74083,1.74083,0,0,1,.62988.46386Z"/>
        </g>'''
    )

    # JAVSCRIPT
    js_content = (pathlib.Path(__file__).parent / "hover_icons.js").read_text(encoding="utf-8")
    js_contents_id = "svg-" + str(uuid.uuid4())
    js_content = js_content.replace("__REPLACE_ME__", js_contents_id) # unique id

    # X
    x = get_layers("X", "cls-8", "cls-9", n_layers_x, x_base, y_base, margin_layer, size_cols_x, size_rows_x, f'({adata.X.shape[0]}, {adata.X.shape[1]}, {n_layers_x})')


    # OBS
    obs = (
        f'<g id="obs" class="block">'
            f'<rect class="cls-2" x={x_base_obs} y={y_base} width={size_cols_obs} height={size_rows_x} />'
            f'<text class="cls-3" x={x_base_obs + size_cols_obs/2} y={y_base + size_rows_x/2}>obs</text>'
            f'<text class="cls-3" x={x_base_obs + size_cols_obs/2} y={y_base + size_rows_x/2+10}>{(adata.obs.shape[0], adata.obs.shape[1])}</text>'
        f'</g>'
    )


    # VAR
    var = (
        f'<g id="var" class="block">'
            f'<rect class="cls-11" x={x_base} y={y_base_var} width={size_cols_x} height={size_rows_var} />'
            f'<text class="cls-3" x={x_base + size_cols_x/2} y={y_base_var + size_rows_var/2}>var</text>'
            f'<text class="cls-3" x={x_base + size_cols_x/2} y={y_base_var + size_rows_var/2+10}>{(adata.var.shape[0], adata.var.shape[1])}</text>'
        f'</g>'
    )


    # OBSM
    obsm = get_layers("obsm", "cls-7", "cls-7", n_layers_obsm, x_base_obsm, y_base, margin_layer, size_cols_obsm, size_rows_x)


    # OBSP
    obsp = get_layers("obsp", "cls-4", "cls-4", n_layers_obsp, x_base_obsp, y_base, margin_layer, size_rows_x, size_rows_x)


    # VARM 
    varm = get_layers("varm", "cls-6", "cls-6", n_layers_varm, x_base, y_base_varm, margin_layer, size_cols_x, size_rows_varm)


    # VARP
    varp = get_layers("varp", "cls-5", "cls-5", n_layers_varp, x_base, y_base_varp, margin_layer, size_cols_x, size_cols_x)


    # UNS
    uns = (
        '''<g id="Uns">
            <use xlink:href="#uns-braces-left" class="uns-stroke"/>
            <use xlink:href="#uns-braces-left" class="uns-braces"/>
            <use xlink:href="#uns-braces-dots" class="uns-stroke"/>
            <use xlink:href="#uns-braces-dots" class="uns-braces"/>
            <use xlink:href="#uns-braces-right" class="uns-stroke"/>
            <use xlink:href="#uns-braces-right" class="uns-braces"/>
            <use xlink:href="#uns-name-u" class="uns-stroke"/>
            <use xlink:href="#uns-name-u" class="uns-name"/>
            <use xlink:href="#uns-name-n" class="uns-stroke"/>
            <use xlink:href="#uns-name-n" class="uns-name"/>
            <use xlink:href="#uns-name-s" class="uns-stroke"/>
            <use xlink:href="#uns-name-s" class="uns-name"/>
        </g>'''
    )


    # COMBINE
    svg = (
        f'<script type="module">{js_content}</script>'
        f'<svg id={js_contents_id} width={x_total} height={y_total}>'
            f'<defs>'
                f"{style}"
                f"{uns_specification}"
            f'</defs>'
            f'<g id="Layer_1">'
                f'{x}'
                f'{obs}'
                f'{var}'
                f'{obsm}'
                f'{obsp}'
                f'{varm}'
                f'{varp}'
                f'{uns}'
            f'</g>'
        f"</svg>"
    )

    return svg




def _get_ratio_draw(ratio):  
    if ratio > 0.3: 
        ratio_draw = ratio
    elif ratio > 0.05:
        ratio_draw = math.log(1.05 + ratio)
    else: 
        ratio_draw = math.log(1.1)
    return ratio_draw

def get_size_x(n_rows_x, n_cols_x, width, height, ratio_of_canvas = 0.3):
    ratio = n_rows_x / n_cols_x
    reverse = False

    if ratio > 1: 
        reverse = True
        ratio = 1/ratio

    ratio_draw = _get_ratio_draw(ratio)

    if reverse: 
        size_rows_x = ratio_of_canvas * height
        size_cols_x = ratio_draw * ratio_of_canvas * width
    else: 
        size_rows_x = ratio_draw * ratio_of_canvas * height
        size_cols_x = ratio_of_canvas * width

    return size_rows_x, size_cols_x

def get_layers(g_name, class_name_front, class_name_back, n_layers, x_base, y_base, margin_layer, width, height, shape_text=""):
    layer_list = []
    if n_layers > 0:
        if n_layers > 1:
            for i in range(n_layers-1,0,-1): 
                layer_list.append(f'<rect class={class_name_back} x={x_base - i*margin_layer} y={y_base + i*margin_layer} width={width} height={height} />')
        layer_list.append(f'<rect class={class_name_front} x={x_base} y={y_base} width={width} height={height} />')
    layers = (
        f'<g id={g_name} class="block">'
            f'{"".join(layer_list)}'
            f'<text class="cls-3" x={x_base + width/2} y={y_base + height/2}>{g_name}</text>'
            f'<text class="cls-3" x={x_base + width/2} y={y_base + height/2+10}>{shape_text}</text>'
        f'</g>'
    )
    return layers
