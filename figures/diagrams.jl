# Self-contained SVG redraws of the method illustrations cited by plasma.tex.
# The drawings explain the algorithms; they are not numerical data or EPS copies.
# Run from the repository root: julia --project=. figures/diagrams.jl

using Printf

const OUT = joinpath(@__DIR__, "out")
const INK = "#183047"
const BLUE = "#2866a0"
const RED = "#d54c40"
const GREEN = "#21866f"
const GOLD = "#db9a32"
const PALE = "#f3f7fa"
const MUTED = "#5d7180"

escape_xml(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")

function line(io, x1, y1, x2, y2; color = INK, width = 2, dash = "", arrow = false)
    marker = arrow ? " marker-end=\"url(#arrow)\"" : ""
    style = isempty(dash) ? "" : " stroke-dasharray=\"$dash\""
    print(io, "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" stroke=\"$color\" stroke-width=\"$width\"$style$marker/>\n")
end

function rect(io, x, y, w, h; fill = "white", stroke = "none", width = 1, radius = 12)
    print(io, "<rect x=\"$x\" y=\"$y\" width=\"$w\" height=\"$h\" rx=\"$radius\" fill=\"$fill\" stroke=\"$stroke\" stroke-width=\"$width\"/>\n")
end

function circle(io, x, y, r; fill = BLUE, stroke = "white", width = 2)
    print(io, "<circle cx=\"$x\" cy=\"$y\" r=\"$r\" fill=\"$fill\" stroke=\"$stroke\" stroke-width=\"$width\"/>\n")
end

function label(io, x, y, s; size = 18, color = INK, anchor = "start", weight = 400)
    print(io, "<text x=\"$x\" y=\"$y\" fill=\"$color\" font-family=\"Arial, Helvetica, sans-serif\" font-size=\"$size\" font-weight=\"$weight\" text-anchor=\"$anchor\">$(escape_xml(s))</text>\n")
end

function polyline(io, points; color = BLUE, width = 3, dash = "", fill = "none")
    coords = join((@sprintf("%.1f,%.1f", x, y) for (x, y) in points), " ")
    style = isempty(dash) ? "" : " stroke-dasharray=\"$dash\""
    print(io, "<polyline points=\"$coords\" fill=\"$fill\" stroke=\"$color\" stroke-width=\"$width\" stroke-linecap=\"round\" stroke-linejoin=\"round\"$style/>\n")
end

function particle(io, x, y, charge; radius = 16)
    color = charge == "+" ? RED : BLUE
    circle(io, x, y, radius; fill = color)
    label(io, x, y + 6, charge; size = 21, color = "white", anchor = "middle", weight = 700)
end

function panel(io, x, title; width = 268)
    rect(io, x, 132, width, 305; fill = "white", stroke = "#dce5ec")
    label(io, x + 18, 168, title; size = 18, weight = 700)
end

function draw(content, name, title, subtitle, note)
    mkpath(OUT)
    io = IOBuffer()
    print(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 960 540\" role=\"img\" aria-label=\"$(escape_xml(title))\">\n")
    print(io, "<defs><marker id=\"arrow\" markerWidth=\"10\" markerHeight=\"8\" refX=\"9\" refY=\"4\" orient=\"auto\"><path d=\"M 0 0 L 10 4 L 0 8 z\" fill=\"$INK\"/></marker></defs>\n")
    rect(io, 0, 0, 960, 540; fill = PALE, radius = 0)
    label(io, 48, 51, title; size = 29, weight = 700)
    label(io, 48, 81, subtitle; size = 16, color = MUTED)
    content(io)
    label(io, 48, 500, note; size = 15, color = MUTED)
    print(io, "</svg>\n")
    path = joinpath(OUT, name * ".svg")
    write(path, String(take!(io)))
    return path
end

function tiled_box(io, x, y, size; strong = false, charges = true)
    rect(io, x, y, size, size; fill = strong ? "#e9f3fb" : "white",
         stroke = strong ? BLUE : "#bdccd7", width = strong ? 3 : 1, radius = 0)
    charges || return
    for (u, v, c) in ((0.23, 0.3, "+"), (0.68, 0.18, "−"),
                       (0.75, 0.72, "+"), (0.32, 0.73, "−"))
        particle(io, x + u * size, y + v * size, c; radius = size < 110 ? 10 : 13)
    end
end

function draw_intro()
    draw("intro", "The simulated plasma", "The three ingredients of each time step",
         "Red: fixed ions   ·   Blue: moving electrons   ·   Arrow: laser field") do io
        rect(io, 75, 135, 445, 300; fill = "white", stroke = BLUE, width = 3)
        for (x, y) in ((145, 195), (302, 172), (415, 313), (205, 362))
            particle(io, x, y, "+")
        end
        for (x, y) in ((195, 250), (350, 245), (440, 198), (288, 355))
            particle(io, x, y, "−")
        end
        line(io, 590, 225, 835, 225; color = GOLD, width = 8, arrow = true)
        label(io, 710, 196, "Eₓ(t)", size = 24, color = GOLD, anchor = "middle", weight = 700)
        label(io, 592, 326, "Periodic box", size = 22, weight = 700)
        label(io, 592, 361, "positions + momenta", size = 18)
    end
end

function draw_field()
    draw("field", "The laser is switched on gradually", "Relaxation, linear ramp, then full oscillation",
         "The field is uniform in space and points along x.") do io
        rect(io, 80, 128, 800, 320; fill = "white", stroke = "#dce5ec")
        for (x, w, c) in ((105, 220, "#edf1f4"), (325, 220, "#fff3db"), (545, 300, "#eaf5ef"))
            rect(io, x, 208, w, 172; fill = c, radius = 0)
        end
        line(io, 105, 294, 845, 294; color = MUTED, width = 1)
        line(io, 105, 394, 845, 394; color = INK, width = 2, arrow = true)
        points = [(x, 294 - 82 * clamp((x - 325) / 220, 0, 1) * cos((x - 325) * 2π / 115)) for x in 105:4:840]
        polyline(io, points; color = BLUE, width = 3)
        for (x, s) in ((215, "relax"), (435, "ramp"), (695, "full field"))
            label(io, x, 170, s; size = 21, anchor = "middle", weight = 700)
        end
        label(io, 105, 427, "0", size = 16)
        label(io, 835, 427, "time", size = 16, anchor = "end")
    end
end

function draw_periodic(stage)
    titles = ("One simulation box", "Periodic copies", "Minimum-image separation")
    subtitles = ("Electrons can cross the boundary", "Every charge has infinitely many images",
                 "The nearest image gives the pair displacement")
    draw("periodic_$stage", titles[stage], subtitles[stage],
         "The Ewald correction accounts for the other periodic images.") do io
        if stage == 1
            tiled_box(io, 260, 145, 320; strong = true)
            line(io, 480, 317, 650, 317; color = BLUE, width = 4, arrow = true)
            label(io, 658, 324, "wrap", size = 20, color = BLUE)
        else
            for row in 0:2, col in 0:2
                tiled_box(io, 300 + col * 112, 112 + row * 112, 112;
                          strong = row == 1 && col == 1, charges = stage == 2)
            end
            if stage == 3
                circle(io, 488, 300, 10; fill = BLUE)
                circle(io, 550, 300, 10; fill = RED)
                line(io, 501, 300, 536, 300; color = GREEN, width = 4, arrow = true)
                label(io, 660, 294, "short periodic", size = 18, color = GREEN)
                label(io, 660, 319, "displacement", size = 18, color = GREEN)
            end
        end
    end
end

function draw_ewald(stage)
    if stage == 1
        draw("ewald_1", "Ewald splitting: a periodic charge", "Add and subtract a smooth Gaussian around every point charge",
             "The local part has zero net charge; the smooth part has a Fourier sum.") do io
            for (x, t) in ((45, "periodic points"), (345, "point − Gaussian"),
                           (645, "Gaussian lattice"))
                panel(io, x, t)
            end
            for (x0, y) in ((65, 360), (365, 360), (665, 360))
                line(io, x0, y, x0 + 230, y; color = MUTED, width = 1)
            end
            for x in 95:50:285
                line(io, x, 355, x, 240; color = RED, width = 4)
            end
            for x in 395:50:585
                points = [(u, 355 + 30 * exp(-((u - x) / 15)^2)) for u in (x - 22):2:(x + 22)]
                polyline(io, points; color = RED, width = 2)
                line(io, x, 355, x, 240; color = RED, width = 3)
            end
            points = [(u, 355 - 65 * sum(exp(-((u - x) / 18)^2) for x in 695:50:885)) for u in 665:3:890]
            polyline(io, points; color = GREEN, width = 3)
            label(io, 321, 303, "=", size = 33, anchor = "middle", weight = 700)
            label(io, 621, 303, "+", size = 33, anchor = "middle", weight = 700)
            label(io, 179, 407, "point-charge comb", size = 16, anchor = "middle")
            label(io, 479, 407, "local, neutral", size = 16, anchor = "middle")
            label(io, 779, 407, "smooth, periodic", size = 16, anchor = "middle")
        end
    else
        draw("ewald_2", "The potential used by the simulation", "The smooth periodic correction is precomputed once",
             "Pair softening is added separately when forces are evaluated.") do io
            for (x, t) in ((45, "periodic potential"), (345, "central 1/r"),
                           (645, "Ewald correction"))
                panel(io, x, t)
            end
            for x in (65, 365, 665)
                line(io, x, 352, x + 230, 352; color = MUTED, width = 1)
            end
            xvalues = 1:3:225
            pair(u) = 343 - 120 / (abs(u - 112) / 12 + 1)
            correction(u) = 312 + 0.002 * (u - 112)^2
            left = [(65 + u, pair(u) + correction(u) - 352) for u in xvalues]
            middle = [(365 + u, pair(u)) for u in xvalues]
            right = [(665 + u, correction(u)) for u in xvalues]
            polyline(io, left; color = BLUE)
            polyline(io, middle; color = RED)
            polyline(io, right; color = GREEN)
            label(io, 321, 292, "=", size = 33, anchor = "middle", weight = 700)
            label(io, 621, 292, "+", size = 33, anchor = "middle", weight = 700)
            label(io, 480, 401, "singular only at r = 0", size = 15, anchor = "middle")
            label(io, 780, 401, "smooth across the box", size = 15, anchor = "middle")
        end
    end
end

const TREE_POINTS = ((112, 222), (179, 285), (239, 203), (306, 352), (351, 251),
                     (400, 310), (545, 192), (590, 238), (632, 284), (675, 182),
                     (713, 351), (765, 245), (820, 324))

function draw_tree(stage)
    if stage == 1
        draw("tree_1", "Build a spatial tree", "Subdivide occupied cells, then store their moments",
             "The actual code uses an octree in three dimensions.") do io
            for (x, t) in ((45, "particles"), (345, "four cells"), (645, "refine busy cells"))
                panel(io, x, t)
                rect(io, x + 35, 190, 200, 200; fill = "#f7fafc", stroke = INK, radius = 0)
            end
            for x in (480, 780)
                line(io, x, 190, x, 390; color = MUTED, width = 1)
                line(io, x - 100, 290, x + 100, 290; color = MUTED, width = 1)
            end
            line(io, 830, 190, 830, 290; color = MUTED, width = 1)
            line(io, 780, 240, 880, 240; color = MUTED, width = 1)
            for (u, v) in ((0.14, 0.20), (0.21, 0.74), (0.40, 0.61), (0.69, 0.15),
                           (0.76, 0.29), (0.82, 0.83))
                for x in (80, 380, 680)
                    circle(io, x + 200u, 190 + 200v, 7)
                end
            end
        end
    else
        draw("tree_2", "Barnes–Hut: near or far?", "Open nearby cells; approximate distant cells by multipoles",
             "The historical runs use θ = 0.5 and include quadrupole moments.") do io
            rect(io, 75, 145, 810, 300; fill = "white", stroke = "#dce5ec")
            rect(io, 105, 185, 238, 220; fill = "#f6f9fc", stroke = BLUE, width = 2, radius = 0)
            line(io, 224, 185, 224, 405; color = "#9eb5c5", width = 1)
            line(io, 105, 295, 343, 295; color = "#9eb5c5", width = 1)
            for (x, y) in ((145, 228), (193, 339), (268, 231), (299, 358))
                circle(io, x, y, 8)
            end
            particle(io, 736, 296, "−"; radius = 20)
            circle(io, 224, 295, 18; fill = GOLD)
            line(io, 244, 295, 708, 295; color = GREEN, width = 3, dash = "8 7")
            label(io, 480, 275, "distance d", size = 20, color = GREEN, anchor = "middle")
            label(io, 224, 423, "cell size s", size = 18, color = BLUE, anchor = "middle")
            label(io, 480, 358, "accept cell when s / d < θ", size = 22, anchor = "middle", weight = 700)
        end
    end
end

function trajectory_points(xoffset = 0, yoffset = 0)
    points = Tuple{Float64,Float64}[]
    for u in 0:0.025:1
        x = 150 + 630u + xoffset
        y = 250 - 105 * sin(2π * u) * exp(-((u - 0.48) / 0.4)^2) + yoffset
        push!(points, (x, y))
    end
    return points
end

function draw_trajectory(stage)
    titles = ("A close encounter bends the path", "A coarse step misses the bend",
              "Refine the difficult interval", "Keep coarse steps elsewhere")
    subtitles = ("The ion changes one electron's motion strongly",
                 "Local error grows near the ion", "Substeps reveal the turn",
                 "Individual steps follow the local difficulty")
    draw("traject_$stage", titles[stage], subtitles[stage],
         "Blue points are computed positions; dashed stretches show coarse estimates.") do io
        rect(io, 85, 135, 790, 310; fill = "white", stroke = "#dce5ec")
        particle(io, 484, 308, "+"; radius = 24)
        pts = trajectory_points()
        polyline(io, pts; color = "#b6c3cd", width = 3)
        idx = stage == 1 ? collect(1:8:length(pts)) :
              stage == 2 ? [1, 9, 33, 41] :
              stage == 3 ? [1, 9, 15, 19, 22, 26, 33, 41] :
                           [1, 9, 15, 18, 20, 22, 26, 33, 41]
        polyline(io, pts[idx]; color = BLUE, width = 3, dash = stage == 2 ? "8 7" : "")
        for i in idx
            circle(io, pts[i][1], pts[i][2], 7; fill = BLUE)
        end
        if stage >= 2
            rect(io, 382, 186, 215, 210; fill = "none", stroke = GOLD, width = 2)
            label(io, 489, 418, "close encounter", size = 17, color = GOLD,
                  anchor = "middle")
        end
    end
end

function draw_trajectory_step(stage)
    titles = ("Compare two estimates", "Split a rejected step", "Freeze stable particles")
    subtitles = ("One RK4 step versus two half-steps", "Only the difficult particle recurses",
                 "Interpolate the quiet path while the close encounter is refined")
    draw("traject_step$stage", titles[stage], subtitles[stage],
         "This is a schematic of the recursive individual-step integrator.") do io
        rect(io, 70, 126, 820, 324; fill = "white", stroke = "#dce5ec")
        particle(io, 478, 330, "+"; radius = 22)
        near = [(120, 220), (270, 263), (385, 276), (455, 265), (510, 240), (570, 198), (795, 175)]
        far = [(120, 340), (345, 350), (570, 357), (795, 362)]
        polyline(io, near; color = BLUE, width = 3)
        polyline(io, far; color = GREEN, width = 3)
        polyline(io, [near[1], near[end]]; color = RED, width = 2, dash = "8 7")
        for p in (near[1], near[end], far[1], far[end])
            circle(io, p[1], p[2], 7; fill = BLUE)
        end
        if stage >= 2
            for p in near[2:end-1]
                circle(io, p[1], p[2], 6; fill = BLUE)
            end
            label(io, 470, 162, "error > ε → split", size = 20, color = RED,
                  anchor = "middle", weight = 700)
        else
            label(io, 470, 162, "compare end positions", size = 20, color = RED,
                  anchor = "middle", weight = 700)
        end
        if stage == 3
            polyline(io, far; color = GREEN, width = 3, dash = "7 6")
            label(io, 750, 405, "stable → interpolate", size = 18, color = GREEN,
                  anchor = "end", weight = 700)
        else
            label(io, 750, 405, "quiet particle", size = 18, color = GREEN, anchor = "end")
        end
    end
end

function main()
    draw_intro()
    draw_field()
    foreach(draw_periodic, 1:3)
    foreach(draw_ewald, 1:2)
    foreach(draw_tree, 1:2)
    foreach(draw_trajectory, 1:4)
    foreach(draw_trajectory_step, 1:3)
    return nothing
end

main()
