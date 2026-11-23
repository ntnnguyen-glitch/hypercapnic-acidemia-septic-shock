# Baseline table, figure data and Figures 1-4.
source("R/00_setup.R")
R   <- readRDS(file.path(OUT_DIR, "results.rds"))
A   <- as.data.table(readRDS(file.path(WORK_DIR, "analysis.rds")))
RAW <- as.data.table(readRDS(file.path(WORK_DIR, "analytic.rds")))
fl  <- jsonlite::read_json(file.path(OUT_DIR, "flow.json"))
GL  <- c(n = "Normocapnia", lo = "Hypocapnia", chc = "Compensated hypercapnia", hca = "Hypercapnic acidemia")
f2  <- function(x, d = 2) formatC(x, format = "f", digits = d)
pv  <- function(p) ifelse(p < 0.001, "<0.001", f2(p, 3))
nn  <- function(z) prettyNum(z, big.mark = ",")
FONT <- "Arial"
th <- theme_classic(base_size = 9, base_family = FONT) +
  theme(axis.text = element_text(colour = "black"), plot.title = element_text(face = "bold", size = 10))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".tiff")), p, width = w, height = h, dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}

## Table 1: baseline characteristics by blood gas category (patients without COPD)
cont <- c(age = "Age, years", paco2 = "PaCO2, mmHg", ph = "pH", hco3 = "Bicarbonate, mmol/L", pf = "PaO2/FiO2, mmHg",
          lact = "Lactate, mmol/L", bun = "Blood urea nitrogen, mg/dL", creat = "Creatinine, mg/dL", na = "Sodium, mmol/L",
          hb = "Hemoglobin, g/dL", wbc = "White blood cells, x10^9/L", plt = "Platelets, x10^9/L", hr = "Heart rate, beats/min",
          map = "Mean arterial pressure, mmHg", rr = "Respiratory rate, breaths/min", sev = "Severity score")
bin <- c(male = "Male sex", imv = "Invasive ventilation at blood gas", niv = "Noninvasive ventilation, first 24 h", nv2 = "Two or more vasopressors, first 24 h",
         cad = "Coronary artery disease", af = "Atrial fibrillation", ckd = "Chronic kidney disease", dm = "Diabetes mellitus",
         liver = "Chronic liver disease", cancer = "Malignancy", pna = "Pneumonia", osa = "Sleep apnea or obesity hypoventilation", hd = "In-hospital death")
tab1 <- function(d) {
  d <- copy(d); d[, `:=`(d28 = as.integer(!is.na(death_h) & death_h <= 672), nv2 = as.integer(!is.na(nvaso) & nvaso >= 2))]
  lv <- c("lo", "n", "chc", "hca")
  out <- list(data.table(Variable = "Patients, n", g1 = nn(sum(d$q == "lo")), g2 = nn(sum(d$q == "n")), g3 = nn(sum(d$q == "chc")), g4 = nn(sum(d$q == "hca")), P = ""))
  for (v in names(cont)) { x <- d[[v]]; if (all(is.na(x))) next
    dd <- if (v == "ph") 2 else 1
    s <- sapply(lv, function(l) { z <- x[d$q == l]; sprintf("%s (%s-%s)", f2(median(z, na.rm = TRUE), dd), f2(quantile(z, .25, na.rm = TRUE), dd), f2(quantile(z, .75, na.rm = TRUE), dd)) })
    lab <- if (v == "sev") ifelse(d$db[1] == "MIMIC-IV", "SOFA score (non-respiratory, non-neurological)", "APACHE IV score") else cont[v]
    out[[length(out) + 1]] <- data.table(Variable = lab, g1 = s[1], g2 = s[2], g3 = s[3], g4 = s[4], P = pv(kruskal.test(x ~ d$q)$p.value)) }
  for (v in c(names(bin), "d28")) { x <- d[[v]]
    s <- sapply(lv, function(l) { z <- x[d$q == l]; sprintf("%d (%s)", sum(z == 1, na.rm = TRUE), f2(100 * mean(z == 1, na.rm = TRUE), 1)) })
    lab <- if (v == "d28") "28-day death" else bin[v]
    out[[length(out) + 1]] <- data.table(Variable = lab, g1 = s[1], g2 = s[2], g3 = s[3], g4 = s[4], P = pv(suppressWarnings(chisq.test(table(x, d$q)))$p.value)) }
  r <- rbindlist(out); setnames(r, c("Variable", "Hypocapnia", "Normocapnia", "Compensated hypercapnia", "Hypercapnic acidemia", "P")); r }
T1 <- list(`MIMIC-IV` = tab1(RAW[copd == 0]))
for (x in names(T1)) fwrite(T1[[x]], file.path(OUT_DIR, paste0("table1_", x, ".csv")))
miss <- RAW[copd == 0, lapply(.SD, function(z) round(100 * mean(is.na(z)), 1)), by = db, .SDcols = c("pf", "hco3", "lact", "alb", "bun", "creat", "hb", "hr", "map", "sev")]
fwrite(miss, file.path(OUT_DIR, "missing_values.csv"))

## Figure 1: flow chart (arrow heads stop short of the boxes)
bx <- function(x, y, lab, w = 3.5, h = 0.62) data.table(x = x, y = y, lab = lab, w = w, h = h)
m <- R[["MIMIC-IV"]]$t1
ncopd <- RAW[copd == 1, .N]
grp <- function(t) sprintf("Hypocapnia: n = %s\nNormocapnia: n = %s\nCompensated hypercapnia: n = %s\nHypercapnic acidemia: n = %s", nn(t$n[2]), nn(t$n[1]), nn(t$n[3]), nn(t$n[4]))
B <- rbind(
  bx(3, 10.4, sprintf("MIMIC-IV v3.1\nFirst ICU stay, age \u226518 years\nn = %s", nn(fl$m0))),
  bx(3, 9.0, sprintf("Sepsis (ICD-9-CM or ICD-10-CM)\nn = %s", nn(fl$m_sep))),
  bx(3, 7.6, sprintf("ICU stay \u226524 h\nn = %s", nn(fl$m_los))),
  bx(3, 6.2, sprintf("Vasopressor from -6 h to +24 h (septic shock)\nn = %s", nn(fl$m_vaso))),
  bx(3, 4.8, sprintf("Arterial blood gas from -6 h to +24 h\nn = %s", nn(fl$m_abg))),
  bx(3, 3.4, sprintf("Without COPD (primary cohort)\nn = %s", nn(R[["MIMIC-IV"]]$n))),
  bx(6.9, 4.1, sprintf("With COPD, n = %s\n(sensitivity and effect\nmodification analyses)", nn(ncopd)), w = 2.3, h = 0.95),
  bx(3, 1.75, grp(m), h = 1.2))
gap <- 0.08
b0 <- B[x == 3][order(-y)]
ar <- data.table(x = 3, xend = 3, y = b0$y[-nrow(b0)] - b0$h[-nrow(b0)] / 2, yend = b0$y[-1] + b0$h[-1] / 2 + gap)
cb <- B[x == 6.9]; side <- data.table(x = 3, y = 4.1, xend = cb$x - cb$w / 2 - gap, yend = 4.1)
p1 <- ggplot() + geom_rect(data = B, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2), fill = "white", colour = "black", linewidth = 0.4) +
  geom_text(data = B, aes(x, y, label = lab), size = 2.5, family = FONT, lineheight = 0.95) +
  geom_segment(data = ar, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.4, arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
  geom_segment(data = side, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.4, arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
  coord_cartesian(xlim = c(1.0, 8.3), ylim = c(1.0, 10.9), expand = FALSE) + theme_void(base_family = FONT)
save_fig(p1, "Figure1", 5.4, 7.2)

## Figure 2: Kaplan-Meier curves with numbers at risk
cols <- c("Normocapnia" = "black", "Hypocapnia" = "#3182bd", "Compensated hypercapnia" = "#31a354", "Hypercapnic acidemia" = "#de2d26")
kmp <- function(x, tag) {
  d <- A[db == x & copd == 0]; d[, g := factor(GL[as.character(q)], levels = names(cols))]
  s <- survfit(Surv(tt, ev) ~ g, d); k <- data.table(t = s$time, surv = s$surv, g = rep(sub("g=", "", names(s$strata)), s$strata))
  k <- rbindlist(lapply(names(cols), function(z) rbind(data.table(t = 0, surv = 1, g = z), k[g == z]))); k[, g := factor(g, levels = names(cols))]
  pl <- ggplot(k, aes(t, 100 * surv, colour = g)) + geom_step(linewidth = 0.6) + scale_colour_manual(values = cols, name = NULL) +
    scale_x_continuous(breaks = seq(0, 28, 7), limits = c(-0.8, 28.8), expand = c(0, 0)) + coord_cartesian(ylim = c(40, 100)) + guides(colour = guide_legend(ncol = 1)) +
    annotate("text", x = 1, y = 45, hjust = 0, size = 2.8, family = FONT, label = sprintf("Log-rank P %s", ifelse(R[[x]]$km < 0.001, "< 0.001", paste("=", f2(R[[x]]$km, 3))))) +
    labs(x = NULL, y = "Survival probability (%)") + th + theme(legend.position = "inside", legend.position.inside = c(0.78, 0.86), legend.key.height = unit(3.2, "mm"))
  ss <- summary(s, times = seq(0, 28, 7), extend = TRUE); rk <- data.table(t = ss$time, n = ss$n.risk, g = factor(sub("g=", "", ss$strata), levels = rev(names(cols))))
  rt <- ggplot(rk, aes(t, g, label = n)) + geom_text(size = 2.4, family = FONT) + scale_x_continuous(breaks = seq(0, 28, 7), limits = c(-0.8, 28.8), expand = c(0, 0)) +
    labs(x = "Days after the index arterial blood gas", y = NULL, title = "Number at risk") + th +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(face = "plain", size = 8.5))
  pl / rt + plot_layout(heights = c(3.2, 1.2)) }
p2 <- kmp("MIMIC-IV", "")
save_fig(p2, "Figure2", 5.0, 5.2)

## Figure 3: PaCO2 dose-response by acidemia status (restricted cubic splines)
d <- A[copd == 0]; M3x <- "age+male+cad+af+ckd+dm+liver+cancer+pna+osa+sevz+imv+niv+nvaso+lpf+lact+bun+creat+na+hb+wbc+plt+hr+map"
fr <- glm(as.formula(paste("ev~ns(paco2,knots=c(35,45,60),Boundary.knots=c(20,100))*acid+", M3x)), binomial, d)
nd <- d[rep(1, 400)]; for (v in names(nd)) if (is.numeric(nd[[v]])) set(nd, j = v, value = median(d[[v]], na.rm = TRUE))
nd$paco2 <- rep(seq(25, 90, length.out = 200), 2); nd$acid <- rep(0:1, each = 200)
pr <- predict(fr, nd, se.fit = TRUE); r0 <- predict(fr, transform(nd[1], paco2 = 40, acid = 0))
sp <- data.table(x = nd$paco2, acid = nd$acid, or = exp(pr$fit - r0), lo = exp(pr$fit - 1.96 * pr$se.fit - r0), hi = exp(pr$fit + 1.96 * pr$se.fit - r0))
rg <- d[, .(a = quantile(paco2, .025), b = quantile(paco2, .975)), by = acid]; sp <- sp[x >= rg$a[match(acid, rg$acid)] & x <= rg$b[match(acid, rg$acid)]]
sp[, grp := factor(acid, 0:1, c("pH \u22657.35", "pH <7.35"))]
p3 <- ggplot(sp, aes(x, or, colour = grp, fill = grp)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) + geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40", linewidth = 0.4) + geom_vline(xintercept = c(35, 45), linetype = 3, colour = "grey40", linewidth = 0.4) +
  scale_colour_manual(values = c("#3182bd", "#de2d26"), name = NULL) + scale_fill_manual(values = c("#3182bd", "#de2d26"), name = NULL) +
  scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) + scale_x_continuous(breaks = seq(30, 90, 10)) +
  labs(x = "PaCO2 on the index arterial blood gas (mmHg)", y = "Adjusted odds ratio for 28-day death") + th + theme(legend.position = "inside", legend.position.inside = c(0.8, 0.15))
save_fig(p3, "Figure3", 5.2, 4.2)

## Figure 4: subgroup forest plot (hypercapnic acidemia vs normocapnia)
sg <- R[["MIMIC-IV"]]$sg
fd <- rbindlist(lapply(names(sg), function(s) { z <- sg[[s]]; lbl <- s
  data.table(mod = lbl, lvl = c("Yes", "No"), n = c(z$nyes, z$nno), or = c(z$yes[1, 1], z$no[1, 1]), lo = c(z$yes[1, 2], z$no[1, 2]), hi = c(z$yes[1, 3], z$no[1, 3]), pint = z$pint) }))
fd <- fd[is.finite(hi)]; fd[, txt := sprintf("%s (%s-%s)", f2(or), f2(lo), f2(hi))]
mods <- unique(fd$mod); nr <- length(mods) * 3; yy <- nr:1
hdr <- data.table(mod = mods, y = yy[seq(1, nr, 3)], t = sapply(mods, function(z) sprintf("%s (P for interaction = %s)", z, f2(fd[mod == z]$pint[1], 3))))
fd[, y := yy[(match(mod, mods) - 1) * 3 + ifelse(lvl == "Yes", 2, 3)]]
fd[, lab := sprintf("    %s (n = %s)", lvl, nn(n))]
yl4 <- c(0.4, nr + 1.4)
pl <- ggplot() + geom_text(data = fd, aes(0, y, label = lab), hjust = 0, size = 2.5, family = FONT) + geom_text(data = hdr, aes(0, y, label = t), hjust = 0, size = 2.5, family = FONT, fontface = "bold") +
  annotate("text", x = 0, y = nr + 1.1, hjust = 0, size = 2.6, family = FONT, fontface = "bold", label = "Subgroup") + scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = yl4) + theme_void(base_family = FONT)
pm <- ggplot(fd, aes(or, y)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey40", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y", linewidth = 0.4) + geom_point(shape = 15, size = 1.8) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4)) + coord_cartesian(xlim = c(0.4, 4.5)) + scale_y_continuous(limits = yl4) +
  labs(x = "Adjusted OR (95% CI)", y = NULL) + th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
pr4 <- ggplot(fd, aes(0, y)) + geom_text(aes(label = txt), hjust = 0, size = 2.5, family = FONT) +
  annotate("text", x = 0, y = nr + 1.1, hjust = 0, size = 2.6, family = FONT, fontface = "bold", label = "OR (95% CI)") + scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = yl4) + theme_void(base_family = FONT)
p4 <- pl + pm + pr4 + plot_layout(widths = c(1.5, 1.1, 0.7))
save_fig(p4, "Figure4", 7.2, 6.4)

saveRDS(list(T1 = T1, miss = miss), file.path(OUT_DIR, "tables.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info.txt"))
