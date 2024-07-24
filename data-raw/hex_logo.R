library(hexSticker)

logo_image <- fs::path("data-raw", "indicspecies_logo_image.png")
sticker(
  logo_image,
  package = "indicspecies", p_size = 17, p_y = 1.50, p_color = "#000080",
  s_x = 0.97, s_y = .86, s_width = .58,
  filename = fs::path("data-raw", "indicspecies.png"),
  #   url = "emf.creaf.cat", u_size = 6, u_color = "#BFD77A", u_y = .2, u_x = 1.2,
  h_fill = "#E3DEDB", h_color = "#000080"
)
