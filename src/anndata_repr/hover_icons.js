const unique_id = "__REPLACE_ME__"
const svg_id = `svg_${unique_id}`
const output_id = `output_${unique_id}`

let output = document.getElementById(output_id)
console.log('outp',output,output_id)
const svg = output.querySelector("#"+svg_id);
console.log('svg',svg,svg_id)

const blocks = svg.querySelectorAll(".block")


for (let b of blocks) {
    b.addEventListener("mouseover", () => {
        change_opacity(blocks, b, 0.5)
    })
    b.addEventListener("mouseout", () => {
        change_opacity(blocks, 0, 1)
    })
    b.addEventListener("click", () => {
        
        const elements = output.getElementsByClassName('ad-section-summary-in')
        for (let e of elements) {
            if(e.dataset.blockname !== b.getAttribute("id")){
                e.checked = false
            } else {
                e.checked = true
            }
        }

        
    })
}

function change_opacity(blocks, except_block, opacity) {
    for (let b of blocks) {
        if (b !== except_block) {
            b.style.opacity = opacity;
            b.style.transitionDuration = "250ms";
        }
    }
}
