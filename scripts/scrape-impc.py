import scrapy
import re
from scrapy.spiders import CrawlSpider, Rule
from scrapy.linkextractors import LinkExtractor


class IMPCSpider(CrawlSpider):
    name = 'impcSpider'
    allowed_domains = ['mousephenotype.org']
    start_urls = ['https://www.mousephenotype.org/impress/pipelines']

    rules = (Rule(LinkExtractor(allow=('impress',), deny=('beta', 'feedback',
                                'login', 'data')), callback='parse_item', follow=True),)

    def parse_item(self, response):
        if re.search(r'parameterontologies', response.url) and \
            response.css('html body div#SiteWrapper div#Content h2 span.procedurekey'
                         '.dark::text').extract() is not None and \
            len(response.css('html body div#SiteWrapper div#Content h2 span.procedurekey'
                             '.dark::text').extract()) > 0:
            yield {
                response.css(
                    'html body div#SiteWrapper div#Content h2 span.procedurekey'
                    '.dark::text').extract()[0]: response.url
            }